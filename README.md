Species Distribution Modeling with BARTs
========================================


README updated 15 Sept 2025.


Project description
-------------------

This repository contains materials (scripts and some data) for a [Physalia short course](https://www.physalia-courses.org/courses-workshops/barts/) demonstrating the uses of Bayesian additive regression tree (BART) methods for species distribution modeling. All scripts are run in R, with specific methods implemented [`dbarts`](https://cran.r-project.org/web/packages/dbarts) with utilities from [`embarcadero`](https://www.github.com/cjcarlson/embarcadero). The first day's material also uses [`gbm3`](https://github.com/gbm-developers/gbm3) and [`brms`](https://doi.org/https://doi.org/10.18637/jss.v080.i01) for demonstration of key concepts. 


Course description
------------------

### Course Overview

The course will introduce and demonstrate the use of Bayesian Additive Regression Tree (BART) methods for species distribution modeling (SDM) and other ecological applications. We will explain how BART modeling works, and identify how BARTs improve over other commonly used SDM methods. Participants will work through all steps necessary for conducting a SDM study with BARTs in R, using utilities in the embarcadero and dbarts packages to select informative environmental predictors, train and evaluate BART models of species occurrence, and use trained models to predict species presence or absence from new data.

### Targeted Audience and Assumed Background

The course is aimed at advanced students, researchers, and practitioners who have some familiarity with species distribution modeling and want to expand their experience and skills. Participants will need a laptop with a webcam and a good internet connection to participate in interactive live sessions. They should be comfortable working in R, particularly the Rstudio environment, and they should be prepared to install specialized packages, edit and write simple scripts, and manage and read in downloaded datasets.

### Learning Outcomes

By the end of the course, participants will
- Understand the structure of BART machine-learning models and how they compare to similar ML methods, especially for species distribution modeling
- Use embarcadero and dbarts utilities to select predictors, and train BART models of species occurrence
- Visualize and interpret predictor effects and interactions in a trained BART model
- Use trained BART SDM models to project species distributions into new regions or times

### Curriculum

- Pre course: Self-guided introduction and installation of necessary packages
- Day 1: Species Distribution Models and Bayesian stats
	- Lecture 1500-1800, 4h practical
	- SDMs overview/review
	- Worked demonstration of SDMs with Boosted Regression Trees using `gbm3` (Topic 1, below)
	- Bayesian stats overview/review
	- Worked demonstration of Bayesian regression with `brms` (Topic 2)
- Day 2: BARTs and embarcadero
	- Lecture 1500-1800, 4h practical
	- BARTs, theory and comparison to other ML methods for SDMs
	- Worked demonstration of an SDM with BARTs, using `embarcadero` --- predictor selection, model training, and prediction with new data (Topic 3)
- Day 3: BART SDM workflow: predictor selection, model evaluation, troubleshooting
	- Lecture 1500-1800, 4h practical
	- Inspecting predictor partial effects and spatial partial effects (Topic 4)
	- Random-intercept BARTs (Topic 5)

Project contents
----------------

- `topics` --- code and slides for worked examples, organized by topic (see below)
- `output` --- materials output by scripts in the worked examples. Saved models amount to more than 1Gb, and are not in the version-controlled repo or posted for download.
- `data` --- data sets needed for the worked examples are more than 100Mb in total, and not posted here, but you can get them from [this Google Drive folder](https://drive.google.com/drive/folders/1JW5iGgF8S2VoYgzjJ5MMXlSP7kWrnto5?usp=share_link).


Module contents
---------------

The code provided should let you dip into any example and take what you want from it. However, this material was developed for a multi-day course, and that experience is best reproduced if you work through it in the same order as presented. In some cases, results of computationally intensive (slow) analyses are saved in an earlier module for use in a later module.

Here are the topics and brief descriptions of their contents.

1. Species distribution modeling with boosted regression trees
	- Preparing presence and pseudoabsence data to train a SDM for Joshua trees, *Yucca brevifolia* and *Y. jaegeriana*
	- Training bosted regression tree (BRT) SDMs from Joshua tree data using `gbm3`
	- Predictor selection
	- Evaluating a trained model
	- Using a trained model to predict from new data

2. Bayesian regression with `brms`

3. Species distribution modeling with BARTs
	- The `varimp.diag()` predictor selection approach
	- Evaluating a trained working model
	- Using a trained model to predict from new data
	- Comparison to results from a BRT

4. Inspecting predictor partial effects
	- Generating and interpreting partial effect plots
	- Generating and interpreting spatial partial effects ("spartials")

5. Random-intercept BART modeling 	
	
