# RNA angles prediction using DSSR


This code is a wrapper code that uses [[DSSR]](#1-4) to get the torsion and pseudo-torsion angles from RNA data.
It computes the torsion angles (`alpha`, `beta`, `gamma`, `delta`, `epsilon`, `zeta`) and the pseudo-torsion angles (`eta`, `theta`, `eta'`, `theta'`, `eta''`, `theta''`) for each nucleotide in the RNA.


## Example 

Here is an example of output using the `csv` format: 
```csv
rank,base,chi,A/S,alpha,beta,gamma,delta,epsilon,zeta,e-z,BI/BII,eta,theta,eta',theta',"eta""","theta""",name
1,A:...1_:[..C]C,-159.5,anti,176.9,45.0,91.9,-150.0,-80.4,-69.7,BI,,NA,-141.3,NA,-141.9,NA,-116.7,rp01
2,A:...2_:[..C]C,-157.6,anti,-56.5,169.0,44.4,74.3,-155.0,-76.2,-78.8,BI,165.0,-154.5,176.5,-150.3,-143.4,-128.8,rp01
3,A:...3_:[..G]G,-155.7,anti,-69.5,177.4,52.5,78.0,-154.9,-75.3,-79.7,BI,168.4,-144.3,177.4,-143.6,-151.9,-108.6,rp01
4,A:...4_:[..C]C,-160.2,anti,-60.3,170.9,42.3,81.8,-154.3,-70.0,-84.4,BI,163.5,-147.7,175.0,-147.7,-134.0,-131.2,rp01
```

## Installation

You can either install the code locally or using Docker. 

If you want to install it locally, you need to have `gcc` installed and python (I used version `3.10.9`). Then, run: 

```bash
make install_all
```

OR the equivalent commands: 
```bash 
make -C dssr/src
pip install -r requirements.txt
```

If you want to install it using docker, use: 
```bash
docker build -t rna_angles_prediction_dssr
```

## Run 

Once the installation is done, you can run the code using either Docker or python.

Locally, the command will be:
```bash
python -m src.dssr_wrapper [--input_path] [--output_path] [--to_csv]
```

If you want to use docker, use:
```bash 
docker run -it rna_angles_prediction_dssr [--input_path] [--output_path] [--to_csv]
```
Please note that you would need to use volumes to mount the input and output directories.
See the `Makefile` that shows an example.

## References 

<a id="1">[1]</a> 
Xiang-Jun Lu & Wilma K. Olson (2003). 
‘3DNA: a software package for the analysis, rebuilding and visualization of three-dimensional nucleic acid structures’, 
Nucleic Acids Res. 31(17), 5108-21. 



<a id="2">[2]</a> 
Xiang-Jun Lu & Wilma K. Olson (2008). 
‘3DNA: a versatile, integrated software system for the analysis, rebuilding and visualization of three-dimensional nucleic-acid structures’, 
Nat. Protoc. 3(7), 1213-27.


<a id="3">[3]</a> 
Guohui Zheng, Xiang-Jun Lu & Wilma K. Olson (2009). 
‘Web 3DNA—a web server for the analysis, reconstruction, and visualization of three-dimensional nucleic-acid structures’, 
Nucleic Acids Res. 37 (Web Server issue), W240–W246. 


<a id="4">[4]</a> 
Andrew Colasanti, Xiang-Jun Lu & Wilma K. Olson (2013). 
‘Analyzing and building nucleic acid structures with 3DNA. Journal of visualized experiments’, 
JoVE, 74, e4401.
