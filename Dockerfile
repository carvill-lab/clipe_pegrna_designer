FROM python:3.11

RUN addgroup --system app && adduser --system --ingroup app app
WORKDIR /home/app
RUN chown app:app -R /home/app

RUN apt-get update && apt-get -y install rsync

COPY src/requirements.txt .
RUN pip install --no-cache-dir --upgrade -r requirements.txt

COPY src/preprocessing preprocessing
RUN python preprocessing/download_genome.py

ADD "https://www.random.org/cgi-bin/randbyte?nbytes=10&format=h" skipcache
RUN python preprocessing/download_clinvar.py

COPY src .

USER app
EXPOSE 8080

CMD ["uvicorn", "app:app", "--host", "0.0.0.0", "--port", "8080"]
