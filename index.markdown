---
# Feel free to add content and custom Front Matter to this file.
# To modify the layout, see https://jekyllrb.com/docs/themes/#overriding-theme-defaults

layout: home
---

# CONS_V 503-A_101 - Topics in Conservation (A) - Fall 2024

## Description 

The purpose of this course is to provide graduate students with the theoretical knowledge and practical skills for the evolutionary analysis of high throughput genomic sequencing data. The course will involve data retrieval and genome assembly, sequence alignment, variant calling and data analysis. The course will be presented as a series of short lectures and lab exercises over a one week period in November 2024.

This year, the workshop runs from November 4th 2024 - November 8th. There are two sessions a day, from 10am to Noon and from 2pm to 4pm.

The course will be run in person on the UBC campus in Biological Sciences building Rm 2139.

# Instructors
Tom Booker
Bernadette Johnson

# Format
A mix of lecture and lab exercises, running in a 2-hour block.

The course material is organized in several topics, with slides and coding examples given in each. Typically we have a lecture for the first half of each class and a practical session for the second.

# Prequisites
Students should have basic knowledge in R and some command line knowledge (although the latter could be obtained during the course)

To get up to speed on working with a Unix system, we recommend having taken the equivalent of Introduction to the Command Line for Genomics from Data Carpentry.

For a handy reference of common Unix commands check out the unix help file. There are some resources there that will help you find the specific command you need for each task.

Evaluation
Your mark will reflect participation in class (discussions and lab exercises etc) as well performance on a final assignment. The breakdown is 80% of your mark will come from participation and 20% from the assignment.

Assignments
At the end of the week you will be given details of an assignment in class that combines different aspects of the work you’ll do throughout the week.

Syllabus
Topic 1 Broad introduction: Scope of course, goals, overview of technology and bioinformatics [Booker,T]
Topic 2 Programming Basics [Booker, T]
Topic 3 Sequence Data and Bioinformatics Gotchyas [Booker, T]
Topic 4 Genome Assembly [Kreiner, J]
Topic 5 RNAseq + differential expression analysis [Booker, T]
Topic 6 Alignment: algorithms and tools [Booker, T]
Topic 7 SNP and variant calling [Kreiner, J]
Topic 8 Patterns of population structure [Kreiner, J]
Topic 9 Signals of adaptation [Kreiner, J]
Topic 10 Data visualization in R [Kreiner, J]
Extras
In previous years, there have been other topics that were included in the workshop as well as brief tutorials on other aspects of bioinformatics. In this section we include the slides and materials for those topics for those who are interested:

Phylogenetic analysis

Chinook
In this workshop we make use of simulated datasets for all of the tutorials and demonstrations. The simutions model a populations of Chinook Salmon living in the Fraser River in Southern British Columbia. We use a population genetic simulation modelling local adaptation to varying environmental conditions. These simulations are used to generate DNA and RNA sequence data that we use to demonstrate bioinformatic principles and give workshop participants experience using standard tools. The biggest benefit to using simulations is that we can compare bioinformatic estimates to a ground truth, which is not possible when using real data.

The GitHub page for Chinook is here: https://github.com/TBooker/Chinook

Obtaining all the files on this site
You may use your internet connection to browse this site, or you may download the entirety of the files on the site in one constantly updated zip archive here

This method dosesn’t require git, however, you’ll have to manually update the files this way (by downloading the whole repo again).

To obtain to all the files via git, type:

git clone https://github.com/UBC-biol525D/UBC-biol525D.github.io.git
To update the all the files at any point in the future, change to the biol525D directory that was created by the previous command and type:

git pull
Use and modification of these resources
You may use any of the materials provided here, and modify them in any way, provided there is appropriate attribution according the license found below and included with this project.

License and Copyright
Copyright (C) 2015 S. Evan Staton, Sariel Hubner, Sam Yeaman

Modified work (c) 2016, 2017, 2018 Gregory Owens, Kathryn Hodgins

Modified work (c) 2019 Gregory Owens, Kathryn Hodgins, J.S. Legare

Modified work (c) 2020 Kathryn Hodgins, Julia Kreiner, Tom Booker

Modified work (c) 2021 Tom Booker, Julia Kreiner

Modified work (c) 2022 Tom Booker, Julia Kreiner

This program is distributed under the MIT (X11) License, which should be distributed with the package. If not, it can be found here: http://www.opensource.org/licenses/mit-license.php

About this site:
This site is powered by GithubPages, and the code backing it is on GitHub here.