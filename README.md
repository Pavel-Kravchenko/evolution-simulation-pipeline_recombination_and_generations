# evolution-simulation-pipeline_recombination_and_generations

<script type="text/javascript">
var s=[],s_timer=[];
 function show(id,h,spd)
 { s[id]= s[id]==spd? -spd : spd;
 s_timer[id]=setTimeout(function() {
 var obj=document.getElementById(id);
 if(obj.offsetHeight+s[id]>=h){obj.style.height=h+"px";obj.style.overflow="auto";}
 else if(obj.offsetHeight+s[id]<=0){obj.style.height=0+"px";obj.style.display="none";}
 else {obj.style.height=(obj.offsetHeight+s[id])+"px";
 obj.style.overflow="hidden1";
 obj.style.display="block";
 setTimeout(arguments.callee, 10);
 }
 }, 10);
 }
</script> 

This is the evolution simulator to compute populations with all possible ranges of recombination, generations and mutations. It creates sets of sequences with different evolution parameters

This pipeline is the part of bioengineering and bioinformatics faculty coursework


## Before you start

The pipeline is available only for <i>Linux</i> users </br>
Make sure that you have installed all components:
<ul>
<li>Python 3.6 or upper https://www.python.org/
<li>Biopython 1.70 or upper http://biopython.org/
<li>R and R-Studio with basic packages https://www.rstudio.com/
</ul>


## Getting started

### Installation

First of all you have to ```clone``` this directory</br></br>
```git clone https://github.com/Pavel-Kravchenko/evolution-simulation-pipeline_recombination_and_generations/```</br></br>
Then ```cd``` in evolution-simulation-pipeline_Ne_and_generations</br></br>
```cd evolution-simulation-pipeline_recombination_and_generations```</br></br>
And ```tar``` Species_files.tar.gz</br></br>
```tar -xvf evolution-simulation-pipeline_recombination_and_generations.tar.gz```</br></br>
Command ```ls -1``` and make sure that you have all files in your directory
```
alignment_reader_script.py
evolution-simulation-pipeline_recombination_and_generations.tar.gz
heatmap_df_maker.py
README.md
R_plot.R
run_script.sh
simulation_script.py
```
Now you are ready to start.
To check, command 
```bash run_script.sh``` and wait until the program is completed. </br></br>
!!!You may have to wait a couple of hours!!!</br></br>
It has to create ``./Files`` directory and archive it to ``archive.tar``
You are expected to receive such demo results in ./Files:

```
out_log_file_1.txt
out_log_file_2.txt
out_log_file_4.txt
simulation_heteroplasmy|gen_128|het-rank_0|mutations_1_0.25
simulation_heteroplasmy|gen_128|het-rank_0|mutations_2_0.25
simulation_heteroplasmy|gen_128|het-rank_0|mutations_4_0.25
simulation_heteroplasmy|gen_128|het-rank_100|mutations_1_0.25
simulation_heteroplasmy|gen_128|het-rank_100|mutations_2_0.25
simulation_heteroplasmy|gen_128|het-rank_100|mutations_4_0.25
simulation_heteroplasmy|gen_128|het-rank_10|mutations_1_0.25
simulation_heteroplasmy|gen_128|het-rank_10|mutations_2_0.25
simulation_heteroplasmy|gen_128|het-rank_10|mutations_4_0.25
simulation_heteroplasmy|gen_128|het-rank_1|mutations_1_0.25
simulation_heteroplasmy|gen_128|het-rank_1|mutations_2_0.25
simulation_heteroplasmy|gen_128|het-rank_1|mutations_4_0.25
simulation_heteroplasmy|gen_128|het-rank_25|mutations_1_0.25
simulation_heteroplasmy|gen_128|het-rank_25|mutations_2_0.25
simulation_heteroplasmy|gen_128|het-rank_25|mutations_4_0.25
simulation_heteroplasmy|gen_128|het-rank_50|mutations_1_0.25
simulation_heteroplasmy|gen_128|het-rank_50|mutations_2_0.25
simulation_heteroplasmy|gen_128|het-rank_50|mutations_4_0.25
simulation_heteroplasmy|gen_128|het-rank_80|mutations_1_0.25
simulation_heteroplasmy|gen_128|het-rank_80|mutations_2_0.25
simulation_heteroplasmy|gen_128|het-rank_80|mutations_4_0.25
simulation_heteroplasmy|gen_16|het-rank_0|mutations_1_0.25
simulation_heteroplasmy|gen_16|het-rank_0|mutations_2_0.25
simulation_heteroplasmy|gen_16|het-rank_0|mutations_4_0.25
simulation_heteroplasmy|gen_16|het-rank_100|mutations_1_0.25
simulation_heteroplasmy|gen_16|het-rank_100|mutations_2_0.25
simulation_heteroplasmy|gen_16|het-rank_100|mutations_4_0.25
simulation_heteroplasmy|gen_16|het-rank_10|mutations_1_0.25
simulation_heteroplasmy|gen_16|het-rank_10|mutations_2_0.25
simulation_heteroplasmy|gen_16|het-rank_10|mutations_4_0.25
simulation_heteroplasmy|gen_16|het-rank_1|mutations_1_0.25
simulation_heteroplasmy|gen_16|het-rank_1|mutations_2_0.25
simulation_heteroplasmy|gen_16|het-rank_1|mutations_4_0.25
simulation_heteroplasmy|gen_16|het-rank_25|mutations_1_0.25
simulation_heteroplasmy|gen_16|het-rank_25|mutations_2_0.25
simulation_heteroplasmy|gen_16|het-rank_25|mutations_4_0.25
simulation_heteroplasmy|gen_16|het-rank_50|mutations_1_0.25
simulation_heteroplasmy|gen_16|het-rank_50|mutations_2_0.25
simulation_heteroplasmy|gen_16|het-rank_50|mutations_4_0.25
simulation_heteroplasmy|gen_16|het-rank_80|mutations_1_0.25
simulation_heteroplasmy|gen_16|het-rank_80|mutations_2_0.25
simulation_heteroplasmy|gen_16|het-rank_80|mutations_4_0.25
...

```


``out_log_file_1.txt`` contains the Kendall's rank correlation tau table

```
Het\Gen	4	8	16	32	64	128	256	512	
0	NA	NA	NA	NA	NA	NA	NA	NA	
1	NA	NA	NA	-0.535	0.224	NA	NA	NA	
10	NA	NA	NA	NA	-0.138	-0.15	-0.138	-0.215	
25	NA	NA	NA	NA	-0.179	-0.228	-0.206	-0.271	
50	NA	NA	NA	NA	-0.354	-0.259	-0.265	-0.298	
80	NA	NA	NA	NA	-0.183	-0.209	-0.206	-0.166	
100	NA	NA	NA	NA	-0.342	-0.352	-0.232	-0.233	
```

To create your populations and analyze them, change necessary parameters in run_script.sh and heatmap_df_maker.py
You may variate following parameters:

<div><a href="#open"  style="decoration: none" onclick="show('hidden1',500,50)"><center> <p>Отчёт по работе программы FastQC</p>
</center>
</a> </div>
<div id=hidden1 style="display:none;height:200px;width:100%;background-color:#f6f1e9">
</div>

## Contact me

Feel free to contact me for any suggestions or critique.

Email: pavel-kravchenk0@yandex.ru 

Site: http://kodomo.fbb.msu.ru/~pavel-kravchenko/index.html 

