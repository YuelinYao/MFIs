# get shiny serves plus tidyverse packages image
FROM rocker/shiny-verse:4.2.3

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    pandoc-citeproc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    python3 \
    python3-dev \
    python3-venv


# copy the app to the image
COPY app.R app.R
COPY renv.lock renv.lock
COPY data data
COPY Produce_devStates.py Produce_devStates.py
COPY app app
COPY www www
COPY renv renv
COPY README.md README.md
COPY .Rprofile .Rprofile

# install R packages required 
# (change it dependeing on the packages you need)
RUN Rscript -e 'install.packages("renv")'
RUN Rscript -e 'renv::consent(provided = TRUE)'
RUN Rscript -e 'renv::restore()'

RUN Rscript -e 'reticulate::virtualenv_create(envname = "example_env_name", python = "python3")'
RUN Rscript -e 'reticulate::virtualenv_install("example_env_name", packages = c("numpy","pandas","scipy","scikit-learn"), ignore_installed=TRUE)'
RUN Rscript -e 'reticulate::use_virtualenv("example_env_name", required = T)'




# select port
EXPOSE 3838

# allow permission

# run app on container start
CMD ["Rscript", "-e", "shiny::runApp('app.R', host = '0.0.0.0', port = 3838)"]






