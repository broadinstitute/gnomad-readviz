DOCKER_IMAGE_NAME=gnomad-readviz
TAG1 = weisburd/${DOCKER_IMAGE_NAME}:v4
TAG2 = weisburd/${DOCKER_IMAGE_NAME}:latest

PWD=`pwd`

all: build push update_scripts

build:
	cd ..
	docker build -t $(TAG1) -t $(TAG2) -f docker/Dockerfile .

push:
	docker push $(TAG1)
	docker push $(TAG2) | tee /dev/tty | grep sha256 | grep latest | cut -d ' ' -f 3 > sha256.txt

update_scripts:
	cat sha256.txt && [ ! -z "`cat sha256.txt`" ] && sed -i.bak "s/${DOCKER_IMAGE_NAME}@sha256:[^\"]*/"${DOCKER_IMAGE_NAME}@`cat sha256.txt`"/"  ../*.py
	rm  ../*.py.bak
