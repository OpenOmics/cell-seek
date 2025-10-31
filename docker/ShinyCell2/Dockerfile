FROM rocker/shiny:4.4.3
LABEL maintainer="Skyler Kuhn <kuhnsa@nih.gov>, Ryan Routsong <routsongrm@nih.gov>"
RUN mkdir -p /opt2/work && mkdir -p /data2 
WORKDIR /opt2
ENV TZ=America/New_York
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime \
        && echo $TZ > /etc/timezone
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update \
    && apt-get upgrade -y \
    && apt-get install -y \
        build-essential \
        bzip2 \
        ca-certificates \
        cmake \
        curl \
        gawk \
        git \
        gzip \
        libfreetype6-dev \
        libgit2-dev \
        libpng-dev \
        libtiff5-dev \
        libjpeg-dev \
        libharfbuzz-dev \
        libfribidi-dev \
        libhdf5-dev \
        libxml2-dev \
        locales \
        make \
        python3 \
        python3-pip \
        unzip \
        wget \
        zlib1g-dev \
        pigz \
        vim \
    && apt-get clean && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
RUN localedef -i en_US -f UTF-8 en_US.UTF-8
RUN ln -sf /usr/bin/python3 /usr/bin/python
RUN mkdir -p /opt2/ \
        && wget --no-check-certificate -O /opt2/which.tar.gz https://ftp.gnu.org/gnu/which/which-2.21.tar.gz \
        && tar -xvf /opt2/which.tar.gz -C /opt2/ \
        && rm /opt2/which.tar.gz \
        && cd /opt2/which-2.21/ \
        && ./configure --prefix=/opt2/which-2.21 \
        && make \
        && make install
COPY install_packages.R /opt2/install_packages.R
RUN cd /opt2 && Rscript /opt2/install_packages.R 
COPY git_installs.sh /opt2/work/git_installs.sh
RUN cd /opt2/work; ./git_installs.sh
RUN Rscript -e 'remotes::install_github("OpenOmics/ShinyCell2@deg_page")'
ADD build_shinycell.R /usr/bin/build_shinycell.R
RUN chmod +x /usr/bin/build_shinycell.R