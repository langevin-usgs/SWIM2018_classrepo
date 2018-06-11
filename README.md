<img src="./notebook/img/header.png" alt="SWIM2018Course" style="width:50;height:20">

# SWIM2018 Pre-Conference Short Course GitHub Repository

This is the class GitHub repository for the pre-conference workshop taught at the [2018 SWIM](https://swim2018.syskonf.pl/).

Pre-Conference Workshop - Modeling Groundwater Flow in Coastal Zones
June 14 – 16, 2018, Gdańsk, Poland.
Course [link](https://swim2018.syskonf.pl/course)

The course is being offered by the conference and will be held at the Mercure Gdańsk Stare Miasto hotel, the venue for SWIM 2016.

Topics will include:
* Theory of variable-density flow and solute transport
* Sharp-interface analytical solutions
* FloPy, a Python tool for the MODFLOW code family
* MODFLOW SWI package
* SEAWAT

## Agenda
[Course Program](https://swim2018.syskonf.pl/conf-data/SWIM2018/files/SWIM%20Short%20Course%20-%20draft%20program.pdf)

The revised agenda is as follows:

### Day 1
* Demonstration -- Introduction to Python and Jupyter Notebooks
* Notebook -- Hand calculations of head and pressure (exHandCalculations_A)
* Presentation -- Analytical solutions and more
* Notebook -- Interface flow toward the coast (exAnalytic_A)
* Notebook -- More interface flow toward the coast (exAnalytic_B)
* Notebook -- A well near the coast (exAnalytic_C)

### Day 2
* Demonstration -- flopy
* Presentation -- Short Introduction to the SWI Package for MODFLOW
* Notebook -- SWI equivalent to analytic example B
* Notebook Exercise -- SWI equivalent to analytic example C (Strack Solution) or A (island)
* Presentation -- Intro to the saltwater intrusion class problem
* Notebook -- SWI equivalent to SEAWAT exB

### Day 3
* Presentation -- SEAWAT concepts
* Presentation -- Overview of Henry Problem
* Notebook -- Henry Problem (exSEAWAT_A)
* Notebook -- Design, run, and calibrate 2D model (exSEAWAT_B)
* Notebook -- Design and run 3D model (exSEAWAT_C)
* Discussion and wrap up

### Optional
* Notebook -- Hand calculations of head and pressure (exHandCalculations_B)
* Demonstration -- Henry analysis (using henry as a function)
* Notebook -- Solute and heat transport (exSEAWAT_D)
* Presentation -- Real world applications of seawater intrusion problems


## Laptop Requirement
Each course participant is expected to arrive with a laptop computer that has the required software installed and tested according to the instructions presented here.  Laptop computers should be running a standard installation of either the Windows or Macintosh operating systems. Users should arrive with privileged account access, sometimes called a “PR account”, in the event that additional software installation is required.  Please coordinate with your IT group prior to arriving to the class.

## Software Requirements
We ask that you install and test the following software prior to showing up for the class.  Installers are located on a public ftp site.  A separate email was sent with instructions about software installation.  Please use the provided installers so that everyone is using the same version.

## Internet Availability
The classroom will have wireless Internet.  Those requiring access to the specific domains will need to establish a VPN connection.

## Exercises/Notebooks
Many of the beginning Python concepts will be taught using the Jupyter Notebook, which runs Python from a web browser.  We will also be running Python scripts using the command line and several other approaches.  

## Python Tutorial

For those with little or no Python experience, we request that you complete a couple of online tutorials.  The first tutorial was developed by Mark Bakker, an instructor for this class.  The tutorial is available at:
 
http://mbakker7.github.io/exploratory_computing_with_python/

The second tutorial is available from Code Academy:

http://www.codecademy.com/en/tracks/python

The Code Academy tutorial includes a variety of topics, but we recommend that students take the following:

* Python syntax
* Strings & Console Output
* Date and Time
* Conditionals and Control Flow
* Functions
* Python Lists and Dictionaries
* Lists and Functions*
* Loops
* Advanced Topics in Python
* Introduction to Classes
* Classes
* File Input/Output


## Repository Folder Structure
The course information is contained in a folder called swim2018shortcourse.  This swim2018shortcourse folder is organized as follows.  This repository (SWIM2018_classrepo.git) will be a folder under the swim2018shortcourse folder.

* installation
* Miniconda3
* software
* SWIM2018_classrepo.git
  * bin
  * doc
  * installation
  * notebook
  * presentation
* working
  
