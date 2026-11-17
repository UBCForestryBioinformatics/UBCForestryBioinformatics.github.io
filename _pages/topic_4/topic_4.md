---
title: "Topic 4 - Genome Annotation"
author: Tom Booker
date: 2024-10-16
category: Jekyll
layout: post
---


# Accompanying material

[Lecture Slides](/pages/topic_4/Topic_4.pdf)


### Install EviAnn from GitHub

You installed FastQC yesterday, but today you're going to install a different program: EviAnn (we discussed it in the lecture)

```
## I lifted this straight from the EviAnn website - you'll need to fix it to work with the current version of the software

tar xvzf EviAnn-X.X.X.tar.gz
cd EviAnn-X.X.X
./install.sh

```



__________
Take a look at the documentation for EviAnn here: [https://github.com/alekseyzimin/EviAnn_release]{https://github.com/alekseyzimin/EviAnn_release}. 



__________ 

*NOTE* EviAnn takes about an 30 mins to run on our servers using 6 threads, so this gives us a good opportunity to run TMUX or BYOBU. Set up EviAnn to run in a TMUX window and leave it running while we discuss what the program does (we'll come back and look at the output later)