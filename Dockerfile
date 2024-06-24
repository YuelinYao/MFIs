# get shiny serves plus tidyverse packages image
FROM rocker/r-base:4.2.3

# system libraries of general use
RUN apt-get update && apt-get install -y \
    sudo \
    pandoc \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    python3 \
    python3-dev \
    python3-venv


# copy the app to the image
COPY . ./


# install R packages required 
RUN Rscript install_Rpackages.R




# select port
EXPOSE 3838

# allow permission

# run app on container start
CMD ["Rscript", "-e", "shiny::runApp('app.R', host = '0.0.0.0', port = 3838)"]






