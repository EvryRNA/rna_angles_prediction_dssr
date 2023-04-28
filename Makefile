build:
	$(MAKE) -C dssr/src
install_python:
	pip install -r requirements.txt

install_all: build install_python

clean:
	$(MAKE) clean -C dssr/src

example_csv:
	python -m src.dssr_wrapper --input_path=data --output_path=rna_puzzles.csv --to_csv

example_json:
	python -m src.dssr_wrapper --input_path=data --output_path=rna_puzzles.json

docker_start:
	docker build -t rna_angles_prediction_dssr .
	docker run -it -v $(PWD)/docker_data:/app/docker_data/ rna_angles_prediction_dssr  --input_path=docker_data/input --output_path=docker_data/output/test.csv --to_csv