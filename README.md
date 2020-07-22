# Analyses for multi round trust game

The model fitting code is based on the code shared by Andreas Hula (https://github.com/AndreasHula).

## How to repduce the analyses

All analyses can be reproduced using the accompanying docker container *mpcoll2/trustgame:latest*
See [docker.io](https://docs.docker.com/get-docker/) on how to install docker on your computer.

After downloading this repository, use your path to this repository to replace $CODEPATH in all commands below.

## Import the data

The *import_prepare_data.py* script loads the raw data, creates long and wide versions of the data frame and saves the data as a binary file for the model fitting. All outputs are saved in the *data* folder.

***TODO : Double check data manipulation in script***

```bash
docker run -it -v $CODEPATH:/code mpcoll2/trustgame:latest python ./code/import_prepare_data.py
```

## Fit the model

The C++ files in the *model_fit* directory are already compiled to run in the docker container. If you are not using the container you might need to modify the scripts and compile them again for your systen. Fitting the model takes 1-2 hours/participant. However, multiple participants can be run in parallel if you have access to more than one cpu thread. RAM usage is about 1GB per participant. To fit the model to one participant or multiple participants, modify the bash script *model_fit_parallel.sh* according to the instructions in the script and run it using the command below.

One output binary file per pair is saved in the *model_fit/outputs* folder.

```bash
docker run -it -v $CODEPATH:/code mpcoll2/trustgame:latest ./code/model_fit_parallel.sh
```

## Read the model parameters

The *read_model_fit.py* script reads the model fitting outputs, and adds the parameters
of the best fitting model for each pair to the raw and wide data files.

```bash
docker run -it -v $CODEPATH:/code mpcoll2/trustgame:latest python ./code/read_model_fit.py
```

## Statistical analyses

**To do **
