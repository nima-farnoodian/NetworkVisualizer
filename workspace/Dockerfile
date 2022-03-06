FROM python:3.11-rc-bullseye

RUN pip3 install virtualenvwrapper
RUN apt update && \
    apt install vim -y

WORKDIR /
COPY ./update-profile.sh ./update-profile.sh
RUN ./update-profile.sh

CMD ["bash"]