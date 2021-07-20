DOCKER_IMAGE_NAME=gnomad-readviz
TAG1 = gcr.io/broad-mpg-gnomad/${DOCKER_IMAGE_NAME}:ukbb

PWD=`pwd`

all: build push update_scripts

build:
	docker build -t $(TAG1) .

push:
	docker push $(TAG1) | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt
