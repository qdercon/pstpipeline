FROM ubuntu:latest

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -qq && apt-get install -y --no-install-recommends  build-essential python3.9 python3-pip python3-setuptools python3-dev wget ca-certificates make g++ software-properties-common dirmngr fonts-open-sans

RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc > /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/"
RUN apt-get update && apt-get install -y --no-install-recommends r-base r-base-dev

WORKDIR /root

COPY requirements.txt /root/requirements.txt
RUN pip3 install -r requirements.txt
	
COPY packages.r /root/packages.r
RUN Rscript /root/packages.r

RUN wget --progress=dot:mega https://github.com/stan-dev/cmdstan/releases/download/v2.28.1/cmdstan-2.28.1.tar.gz && tar -zxpf cmdstan-2.28.1.tar.gz && rm cmdstan-2.28.1.tar.gz
RUN cd cmdstan-2.28.1; make build

COPY notebooks /root/notebooks
COPY data-raw /root/data-raw
COPY outputs /root/outputs

EXPOSE 8888

CMD ["sh", "-c", "jupyter-lab --ip=0.0.0.0 --allow-root"]
