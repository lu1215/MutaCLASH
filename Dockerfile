FROM python:3.5-slim
WORKDIR /MutaCLASH
COPY . .
RUN pip install -r requirements.txt
CMD [ "/bin/bash" ]
