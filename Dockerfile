FROM gcc:9.5.0 as dssr
WORKDIR app
COPY dssr /app/dssr
RUN make -C dssr/src

FROM python:3.10.9
WORKDIR app
COPY --from=dssr /app/dssr /app/dssr
COPY requirements.txt .
RUN pip install -r requirements.txt
# Path to run DSSR
ENV X3DNA=/app/dssr
ENV PATH=/app/dssr/bin:$PATH
COPY . .
ENTRYPOINT ["python", "-m", "src.dssr_wrapper"]