---
title: HW3 - Introduction to R
sidebar: mydoc_sidebar
permalink: mydoc_homework_03.html 
---

## A. Object Subsetting, Import and Export

- __Task 1__: Sort the rows of the `iris` data frame by its first column and sort its columns alphabetically by column names.
- __Task 2__: Subset the first 12 rows, export the result to a text file and view it in a spreadsheet program like Excel or Google Sheets. 
- __Task 3__: Change some column titles in your spreadsheet program, save the result to a tab delimited text file and import it back into R. Note, for this task you only want to include the `read.table` command in the homework result (here R script).

Before you start it can be helpful to evaluate the structure of the `iris` data set with the following commands:
{% highlight r %}
class(iris)
dim(iris)
colnames(iris)
{% endhighlight %}

## B. Scatter Plots
       
- __Task 1__: Generate a scatter plot for the first two columns of the `iris` data frame and color the dots by the `Species` column.
- __Task 2__: Use the `xlim/ylim` arguments to set limits on the x- and y-axes so that all data points are restricted to the bottom left quadrant of the plot. 

Again before you start, evaluate the structure of iris data set. The following commands are useful:

{% highlight r %}
iris[1:4,]
table(iris$Species)
{% endhighlight %}

## C. Bar Plots
        
- __Task 1__: Calculate the mean values for the `Species` components of the first four columns in the `iris` data frame. Organize the results in a matrix where the row names are the unique values from the `iris Species` column and the column names are the names of the first four `iris` columns. 
- __Task 2__: Generate two bar plots for the matrix generated in the previous step: one with stacked bars and one with horizontally arranged bars. 

## D-H. Analysis Worflow

The instructions for these homework assignments are [here](http://girke.bioinformatics.ucr.edu/GEN242/mydoc_Rbasics_15.html).

## Homework submission

Assemble the results from homework assignments A-H in a single R script (`HW3.R`) and upload it to your private GitHub repository under `Homework/HW3/HW3.R`.

## Due date

This homework is due on Thu, April 20th at 6:00 PM.

## Homework Solutions

See [here](https://drive.google.com/open?id=0B-lLYVUOliJFTnRkMFo2SDdHQlE)
