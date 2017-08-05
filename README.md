# Will It Stand
Find out if there is an orientation in which your 3D mesh can stand still

#### Created By:
 - Ziv Shahaf
 - Kujan Lauz

As part of a 3D-Printing Seminar, under the guidance of Professor Gershon Elber, We wanted to identify, given an input 3D model, Which orientations of the model are suitable for 3D printing.
 
To find the orientations we needed to come up with an algorithm for finding all the facets (2D planes) on which the 3D model can "stand still".

## The process

### overview
Start with a 3D mesh input, calculate its convex hul and center of mass.
Using CH and COM find balancing facets


### Mesh Processing
A minimalistic OBJ file parser was created, with some additional implementation of some geometric representation classes: Vector, Matrix, Polygon3D and Geometric tools for 2D and 3D calculations.

<img realwidth="1500" cropx="-36" cropy="49.981" realheight="708" originalsrc="http://lh3.googleusercontent.com/b7WbtYqVTIriuG4ce0r0OHkYJl0_mMmFWm0ULxMJOOqjE2jn-uBxg0pj82j78qrvLtFNd193LVQ47_SoRIEEbzuXDw=s538" src="http://lh3.googleusercontent.com/b7WbtYqVTIriuG4ce0r0OHkYJl0_mMmFWm0ULxMJOOqjE2jn-uBxg0pj82j78qrvLtFNd193LVQ47_SoRIEEbzuXDw=s538" style="position:absolute;left:-36px;top:50px;width:538px;height:254px;">

### Convex Hull Calculation
Used Qhull under copying right (http://www.qhull.org/COPYING.txt) to get the convex hull of the processed mesh. Qhull have a very un-intuitive library/API and lack proper documentation so after finding out how to use it for our needs I answered a [relevant StackOverflow question](https://stackoverflow.com/a/29311240/2523211)

<img realwidth="1500" cropx="-51.5" cropy="42.681" realheight="708" originalsrc="http://lh3.googleusercontent.com/4LxT1gRV9d-ENaJFScmS68XWOzSP1eE03-ScJ_h9YR4C3HAdy4uudoSt9Ea4xK9qlUqK49p6VLGIGo3N4Ep8pmqqddU=s569" src="http://lh3.googleusercontent.com/4LxT1gRV9d-ENaJFScmS68XWOzSP1eE03-ScJ_h9YR4C3HAdy4uudoSt9Ea4xK9qlUqK49p6VLGIGo3N4Ep8pmqqddU=s569" style="position:absolute;left:-51px;top:43px;width:569px;height:269px;">

### Center of Mass Calculation
To calculate the center of mass, we used `volInt.c` (one of its versions can be found [here](https://github.com/OpenFOAM/OpenFOAM-2.1.x/blob/master/src/meshTools/momentOfInertia/volumeIntegration/volInt.c)), the source code of [Fast and Accurate Computation of Polyhedral Mass Properties](http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.56.127&rep=rep1&type=pdf) by Brian Mirtich .


### Finding Balancing Faces
An O(N) algorithm where N is the number of facets of the convex hull.

<img realwidth="1098" cropx="0" cropy="21.968579234972676" realheight="312" originalsrc="http://lh3.googleusercontent.com/GfoXh58w3UrGlpjNLxhUoLq-UhWacGSz3_f0UYKiuUze3Hbl1u3C_7Wi0CXTHe_Zmuagf4bMumWcEUG8XPn2RbZ-=s571" src="http://lh3.googleusercontent.com/GfoXh58w3UrGlpjNLxhUoLq-UhWacGSz3_f0UYKiuUze3Hbl1u3C_7Wi0CXTHe_Zmuagf4bMumWcEUG8XPn2RbZ-=s571" style="position:absolute;left:0px;top:22px;width:571px;height:162px;">

<img realwidth="1500" cropx="-36" cropy="50.012" realheight="708" originalsrc="http://lh3.googleusercontent.com/V9nKBYuLtBMMed5W001cDXhUTePu1arLqP27TrFakVL8mgUw3nIBxrSWyJQaruq4XM2kA3lFF0ojZ7F_lLAKr3ff9o4=s538" src="http://lh3.googleusercontent.com/V9nKBYuLtBMMed5W001cDXhUTePu1arLqP27TrFakVL8mgUw3nIBxrSWyJQaruq4XM2kA3lFF0ojZ7F_lLAKr3ff9o4=s538" style="position:absolute;left:-36px;top:50px;width:538px;height:254px;">
