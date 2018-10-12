# KdTree
 
## Program Structure

The file **parent.py** is provided . It is the program that will calculate the time taken to search. **run.sh** compiles and runs the program.

In order for parent.py to measure the time taken by the program to just answer the query (excluding the k-d tree construction time), the following mechanism has been adopted.

 1. parent.py will execute `sh run.sh <dataset_file>` as a child process. (run.sh should in turn compile and run the program) (Here <**dataset_file**> is the name of the file containing the data points - Format specified below) 
 2. the program (child) constructs the k-d tree using <dataset_file> and output "0" on the standard output (stdout) when done. This would be read by the parent.
 3. parent.py would then write the name/path of the **<query_file>** (which contains the query point for kNN - Format specified below) and the value of **k** (k in kNN) on the standard input (stdin). It also starts the timer now.
 4. the program will now read the name of the <query_file> and the value of k from stdin and process the query using the k-d tree. The output (the k nearest neighbors) should then be written to **results.txt**. (Format specified below)
 5. the program then outputs "1" on stdout so that the parent can stop the timer and check the results.txt.


## File Formats

 - A point is represented as a space separated list of values along each dimension.
 - **<dataset_file>**:
 1. First line is "D N" where D = numer of dimensions and N = number of points. 
 2. N lines follow where each line is a D-dimesnional point.
 - **<query_file>**:
 1. First line is "D" where D = number of dimensions
 2. Next line is the D-dimensional query point.
 - **results.txt**: k lines where each line should be a point.
 
 Example: Suppose the points are 2-dimensional and <dataset_file> contains 3 points, it would look something like
> 2 3  
> 0.0 1.0  
> 1.0 0.0  
> 0.0 0.0

Now if the <query_file> is:
> 2  
> 1.0 1.0

Then for k = 2, results.txt should look like:
> 0.0 1.0  
> 1.0 0.0


## Command to run

    python parent.py <dataset_file> <query_file> <k>
parent.py is compatible with Python 2.x and Python 3.x and with both Linux and Windows. However final evaluation will be done using Python 3.x. So, for students using Python in the program, it is recommended to use Python 3.x.

For any issues, feel free to post on the Piazza group of this course.
