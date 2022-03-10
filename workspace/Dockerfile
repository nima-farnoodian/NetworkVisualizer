FROM python:3.11-rc-bullseye

RUN pip3 install virtualenvwrapper
RUN apt update && \
    apt install vim -y

WORKDIR /nvworkspace
COPY ./update-profile.sh ./update-profile.sh
RUN ./update-profile.sh

EXPOSE 8888

CMD ["bash"]