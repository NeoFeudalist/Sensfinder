{-# LANGUAGE ExistentialQuantification #-}

module Sensfinder where

import Control.Monad
import Control.Monad.Catch
import Data.List
import qualified Numeric.MCMC.NUTS as NUTS
import System.Random
import System.Random.MWC

type Sensitivity = Double
type Mean = Double
type StDev = Double
type LogProb = Double
data Outcome = FirstBetter | SecondBetter
type Datapoint = (Sensitivity, Sensitivity, Outcome)
newtype Dataset = Dataset {points :: [Datapoint]}
data Result = Result {prop1 :: Sensitivity, prop2 :: Sensitivity, optimal :: Sensitivity}

yValue :: Outcome -> Double
yValue FirstBetter = 1
yValue SecondBetter = negate 1

data Config = Config {dataset :: Dataset, meanTheta :: Double, sdTheta :: Double, meanTau :: Double, sdTau :: Double}
data Parameters = Parameters {theta :: Double, tau :: Double} deriving Eq

h :: Double -> Double
h x = if x < 1 then 2**x / 2 else x

f :: Double -> Parameters -> Double
f x params = negate $ (1 / (h $ tau params)) * (x - (h $ theta params)) ^ 2

baseParams :: Parameters
baseParams = Parameters 0 0

opParams :: (Double -> Double -> Double) -> Parameters -> Parameters -> Parameters
opParams op p1 p2 = Parameters (op (theta p1) (theta p2)) (op (tau p1) (tau p2))

addParams :: Parameters -> Parameters -> Parameters
addParams p1 p2 = opParams (+) p1 p2

average :: (Floating a, MonadThrow m) => [a] -> m a
average xs
  | length xs == 0 = throwM $ EmptyListError "xs"
  | otherwise = return $ sum xs / (fromIntegral $ length xs) 

normalLogDensity :: (MonadThrow m) => StDev -> m (Mean -> Double -> LogProb)
normalLogDensity sd
  | sd <= 0 = throwM $ OutOfBoundsErrorGt "sd" 0
  | otherwise = return $ \mean x -> (negate (log sd)) - (log $ 2 * pi) / 2 - sq ((x - mean) / sd) / 2 where
  sq n = n * n

logPrior :: (MonadThrow m) => Config -> m (Parameters -> LogProb)
logPrior conf = do
  nLD1 <- normalLogDensity $ sdTheta conf
  nLD2 <- normalLogDensity $ sdTau conf
  return $ \params -> let
    t1 = nLD1 (meanTheta conf) (theta params)
    t2 = nLD2 (meanTau conf) (tau params)
    in t1 + t2

sigmoid :: (Floating a) => a -> a
sigmoid x = 1 / (1 + exp (- x))

logLikelihood :: Config -> Parameters -> LogProb
logLikelihood conf params = let
  ds = dataset conf
  logPointLikelihood par (a, b, outcome) = let y = yValue outcome in
    log $ sigmoid $ y * ((f a par) - (f b par)) 
      in sum $ map (logPointLikelihood params) (points ds)

logPosterior :: (MonadThrow m) => Config -> m (Parameters -> LogProb)
logPosterior conf = do
  lPrior <- logPrior conf
  let lLikely = logLikelihood conf
  return $ \params -> lPrior params + lLikely params

dNormalLogDensity :: (MonadThrow m) => StDev -> m (Mean -> Double -> Double)
dNormalLogDensity sd
  | sd <= 0 = throwM $ OutOfBoundsErrorGt "sd" 0
  | otherwise = return $ \mean x -> negate $ (x - mean) / (sd * sd)

dLogPrior :: (MonadThrow m) => Config -> m (Parameters -> Parameters)
dLogPrior conf = do
  dnLD1 <- dNormalLogDensity $ sdTheta conf
  dnLD2 <- dNormalLogDensity $ sdTau conf
  return $ \params -> let
    t1 = dnLD1 (meanTheta conf) (theta params)
    t2 = dnLD2 (meanTau conf) (tau params)
      in Parameters t1 t2

dLogLikelihood :: Config -> Parameters -> Parameters
dLogLikelihood conf params = let
  ds = dataset conf
  dLogPointLikelihood par (a, b, outcome) = let
    y = yValue outcome
    hTheta = h $ theta par
    hTau = h $ tau par
    sigmoidTerm = (sigmoid $ negate $ ((a - b) * (a + b - 2 * hTheta) * y) / hTau) - 1
    piecewiseTerm x = if x < 1 then h x * log 2 else 1
    wrtTheta = negate $ (2 * (a - b) * y * sigmoidTerm * piecewiseTerm (theta par)) / hTau
    wrtTau = negate $ ((a - b) * (a + b - 2 * hTheta) * y * sigmoidTerm * piecewiseTerm (tau par)) / (hTau * hTau) in
      Parameters wrtTheta wrtTau 
        in foldr addParams baseParams $ map (dLogPointLikelihood params) (points ds)

dLogPosterior :: (MonadThrow m) => Config -> m (Parameters -> Parameters)
dLogPosterior conf = do
  dlPrior <- dLogPrior conf
  let dlLikely = dLogLikelihood conf
  return $ \params -> addParams (dlPrior params) (dlLikely params)

logPosteriorList :: (MonadThrow m) => Config -> m ([Double] -> LogProb)
logPosteriorList conf = do
  lPosterior <- logPosterior conf
  return $ \params -> lPosterior $ Parameters (params !! 0) (params !! 1)

dLogPosteriorList :: (MonadThrow m) => Config -> m ([Double] -> [Double])
dLogPosteriorList conf = do
  dlPosterior <- dLogPosterior conf
  return $ \params -> let
    res = dlPosterior $ Parameters (params !! 0) (params !! 1)
      in [theta res, tau res]

sample :: Config -> Parameters -> Int -> Int -> IO [Parameters]
sample conf params n drop = do
  when (n < 1) $ throwM $ OutOfBoundsErrorGtEq "n" 1
  when (drop >= n) $ throwM $ DropGtEqSamples drop n
  let chainToParamList = map (\[x0, x1] -> Parameters x0 x1)
  lPosterior <- logPosteriorList conf
  dlPosterior <- dLogPosteriorList conf
  chainToParamList <$> (withSystemRandom . asGenST $
      NUTS.nutsDualAveraging lPosterior dlPosterior n drop [meanTheta conf, meanTau conf])

genInt :: Int -> Int -> IO Int
genInt a b = getStdRandom (randomR (a, b))

sampleList :: [a] -> Int -> IO [a]
sampleList xs n = do
  when (null xs) $ throwM $ EmptyListError "xs"
  when (n <= 0) $ throwM $ OutOfBoundsErrorGtEq "n" 1
  map (xs !!) <$> (replicateM n $ genInt 0 (length xs - 1))

proposal1 :: [Parameters] -> IO Sensitivity
proposal1 = fmap theta . fmap head . flip sampleList 1

proposal2 :: (MonadThrow m) => Sensitivity -> [Parameters] -> m Double
proposal2 p1theta paramChain = average $ map ((subtract p1theta) . (2*) . theta) paramChain

optimalSens :: (MonadThrow m) => [Parameters] -> m Double
optimalSens = average . map theta

runSensfinder :: Config -> Parameters -> Int -> Int -> IO Result
runSensfinder conf params n drop = do
  chain <- sample conf params n drop
  p1 <- proposal1 chain
  p2 <- proposal2 p1 chain
  opt <- optimalSens chain
  if p1 < p2 then
    return $ Result p1 p2 opt
  else
    return $ Result p2 p1 opt

data SensfinderError =
  forall a. (Show a, Num a) => OutOfBoundsErrorGt String a |
  forall a. (Show a, Num a) => OutOfBoundsErrorLt String a |
  forall a. (Show a, Num a) => OutOfBoundsErrorGtEq String a |
  forall a. (Show a, Num a) => OutOfBoundsErrorLtEq String a |
  forall a. (Show a, Num a) => DropGtEqSamples a a |
  EmptyListError String

instance Exception SensfinderError

instance Show Outcome where
  show FirstBetter = "The first sensitivity is better (y = 1)."
  show SecondBetter = "The second sensitivity is better (y = -1)."

instance Show SensfinderError where
  show (OutOfBoundsErrorGt val x) = "Value " ++ val ++ " must be greater than " ++ show x ++ "."
  show (OutOfBoundsErrorLt val x) = "Value " ++ val ++ " must be less than " ++ show x ++ "."
  show (OutOfBoundsErrorGtEq val x) = "Value " ++ val ++ " must be greater than or equal to " ++ show x ++ "."
  show (OutOfBoundsErrorLtEq val x) = "Value " ++ val ++ " must be less than or equal to " ++ show x ++ "."
  show (DropGtEqSamples drop n) = "Cannot drop " ++ show drop ++ " samples out of a total of " ++ show n ++ "."
  show (EmptyListError xsName) = "List " ++ xsName ++ " must not be empty."

instance Show Dataset where
  show ds = let
    pointToLine (a, b, y) = show a ++ ", " ++ show b ++ ", " ++ (show $ yValue y)
      in "Dataset (in this format: sensitivity 1, sensitivity 2, [y = 1.0 if sensitivity 1 is better, y = -1.0 otherwise])" ++ "\n" ++ (intercalate "\n" $ map pointToLine $ points ds)

instance Show Config where
  show conf =
    "Prior theta: Normal(" ++ (show $ meanTheta conf) ++ ", " ++ (show $ sdTheta conf)
    ++ "^2). Prior tau: N(" ++ (show $ meanTau conf) ++ ", " ++ (show $ sdTau conf) 
    ++ "^2)."
    
instance Show Result where
  show res =
    "New proposal: " ++ (show $ prop1 res) ++ " vs. " ++ (show $ prop2 res) ++ ". Optimal sensitivity based on current information: " ++ (show $ optimal res) ++ "."

instance Show Parameters where
  show params = "theta = " ++ (show $ theta params) ++ ", tau = " ++ (show $ tau params) ++ "."

--d = Dataset [(2, 4, 1), (1, 3, negate 1), (2.36, 4.11, 1), (2.44, 3.37, 1), (2.76, 2.88, 1), (1.98, 3.76, 1), (2.82, 3.28, 1)]
d = Dataset [(2,3,FirstBetter), (3,4,SecondBetter)]
--d = Dataset []
c = Config d 2.5 1 1 1
-- c = Config d 1 1 0.2 1 0 10
p = Parameters 0.1 1.2