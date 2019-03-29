# ImageSegmentationGraphCuts
Image segmentation using Graph cuts without network flow overhead
Implements the algorithm proposed by Rachel Silva in her [Thesis](https://scholarworks.rit.edu/theses/9690/)

# To run the program
```
javac *.java
java UI <image_file_name>
```

## UI Menu
* **Scribble**: The options *Draw S* and *Draw T* under the *Scribble* menu can be used to draw the scribbles for source and sink respectively.
* **Compute**: The option *Left Right Most* under *Compute* menu is used to compute the left-most and right-most image segmentations, whereas, the *Prob Sample* option under *Compute* menu is used to find display n probabilistically chosen samples, where n is taken using an input dialog box.
* **Filter**:  The options *Grayscale* and *RGB* under *Filter* can be used to display the image in Gray-scale or RGB mode respectively.

## More details
* [Project Report](https://github.com/rushikv/ImageSegmentationGraphCuts/blob/master/Report.pdf)
* [Rachel Silva's Thesis](https://scholarworks.rit.edu/theses/9690/)
