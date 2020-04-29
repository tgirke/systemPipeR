---
title: News
sidebar: mydoc_sidebar
permalink: mydoc_news.html
---

# *systemPipeR* 1.22 is available

![image](pages/mydoc/miscellaneous_images/SYS_WF.png) 

### OVERVIEW

The following enhancements have been added to systemPipeR.

- With the upgrades provided in this release, systemPipeR has become a much more generic data analysis workflow environment that is no longer limited to analyzing just NGS data. Now it can be efficiently used for data analysis tasks in many omics areas, including genomics, proteomics, metabolomics and drug discovery.

- A workflow control class (`SYSargsList`) has been added allowing users to manage multiple-step workflows from a single container. This way one can select and execute multiple workflow steps with standard R subsetting syntax, e.g. `runWF[1:3]`.

- Various improvements have been added to systemPipeR’s new  command-line interface including the recently introduced `SYSargs2` class that supports the Common Workflow Language ([CWL](https://www.commonwl.org/)). Utilities have been added to visualize workflow designs and topologies with different graphical layouts.

- Improvements have been added to monitor the run status of workflows, as well as tracking of warning and error messages. This includes the generation of both scientific and technical status reports.
<br> <br/>

<script src="https://gist.github.com/dcassol/175fc73f52e62647c45697dbebc81ece.js"></script>