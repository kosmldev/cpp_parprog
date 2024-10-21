FROM gcc:latest
RUN apt-get update && apt-get install -y \
    libopenblas-dev \
    liblapacke-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

COPY ./src /src
WORKDIR /src

RUN g++ main.cpp -lopenblas -llapacke -o program

CMD ["./program"]
