Introduction to shell
================
Haniel Cedraz
2023-05-18

<br>

<br>

First, let’s define what we mean by the “shell.”

In computing, a shell is a program that provides a command-line
interface to interact with the operating system. The shell allows you to
type commands in a terminal and have the computer execute them.

Here are some steps to help you get started with using the shell: Open a
terminal: On Linux, you can use the terminal emulator provided by your
distribution.

Manipulate Files and Directories: You can create, delete, move, and copy
files and directories using commands like mkdir, rm, mv, and cp. For
example, to create a new directory called “my_directory”, you can type
`mkdir my_directory`. To delete a file called “my_file.txt”, you can
type `rm my_file.txt`. Be careful when using these commands, as they can
be powerful and potentially destructive.

First you can use `pwd` to print the current working directory (i.e.,
the directory that the user is currently located in).

``` bash
$ pwd
```

``` bash
$ mkdir my_directory
```

``` bash
$ cp -r my_directory my_directory_test
```

``` bash
$ mv my_directory my_directory_test
```

List Files and Directories: You can use the **`ls`** command to list the
files and directories in the current directory. For example, typing
**`ls`** will list all the files and directories in the current
directory. You can also use the **`-l`** option to display more
information about each file, such as permissions, owner, and size.

``` bash
$ ls my_directory_test   ## Show a more detailed list $ ls -l my_directory_test
```

``` bash
$ cd my_directory
```

Use Pipes and Redirection: You can use pipes and redirection to connect
commands together and manipulate input and output. For example, you can
use the **`|`** symbol to pipe the output of one command as the input to
another command. To redirect the output of a command to a file, you can
use the **`>`** symbol. For example, **`ls > file_list.txt`** will save
the output of the **`ls`** command to a file called “file_list.txt”.

``` bash
$ ls > file_list.txt
```

You can use `touch` to create a new, empty file with a specified name

``` bash
$ touch my_first_file.txt
```

You can use `echo` to display a message or a variable value on the
terminal.

``` bash
$ echo "Hello, World" > my_first_file.txt
```

``` bash
$ echo "I am taking the bioinformatics course" > my_first_file.txt
```

## Text Editor

A text editor is a software application that is used for creating,
editing, and manipulating text files. Text editors are commonly used by
developers, programmers, and technical writers to create and edit code,
scripts, and documentation.

Here are the two most popular text editors in the terminal:

1.  Nano: Nano is a simple and easy-to-use text editor that is ideal for
    beginners. It has a user-friendly interface and provides basic text
    editing features such as syntax highlighting, search and replace,
    and line numbering.

2.  Vim \| Vi: Vim or vi is a powerful and customizable text editor that
    is popular among advanced users and programmers.

``` bash
$ nano my_file_using_nano.txt
```

``` bash
$ vi my_file_using_vi.txt
```

paste this text into the file

``` text
Bioinformatics is a multidisciplinary field that combines biology, computer science, and statistics to analyze and interpret biological data. With the advent of high-throughput technologies, such as next-generation sequencing and microarrays, the amount of biological data generated has exploded, creating a need for powerful computational tools and algorithms to make sense of this data.
```

The **`cat`** command can be used to read files, create new files, and
concatenate files. You can use wildcards to specify groups of files that
match a pattern. For example, you can use **`*`** to match any sequence
of characters, or **`?`** to match a single character. For example,
**`ls *.txt`** will list all files in the current directory that end
with “.txt”.

``` bash
$ cat my_file_using_*
```

Ask ChatGPT to write an assay with 100 itens about bioinformatics.
Create a new file using nano or vi and paste the content in the file.

After saving it you can open the file using the command `less.`

``` bash
$ less your_file.txt
```

Also you can see the first 10 rows or the last 10 rows of the file

``` bash
# The first 10 rows 
$ head your_file.txt  

# You can increase or decrease the number of rows 
$ head -20 your_file.txt 
$ head -60 your_file.txt   

# The last 10 rows 
$ tail your_file.txt  

# You can increase or decrease the number of rows 
$ tail -20 your_file.txt 
$ tail -6 your_file.txt 
```

**`grep`** is used to search for a specific pattern in a file or a set
of files.

``` bash
$ grep "Molecular" your_file.txt
```

``` texinfo
Special Character:
|   Vertical bar, called “pipe”, it indicates that finish a command and start another;
< >     In and out of a command;
>>  Append contente;
.   Current dir
..  Back one dir;
/   Folder Delimiter;
- back to the last dir
```

``` bash
$ cat your_file.txt | grep "Molecular"
```

``` bash
$ cat < your_file.txt | grep "Molecular" > molecular_grep.txt

cat molecular_grep.txt

# Appending in the file
$ cat your_file.txt | grep "genomics" > molecular_grep.txt # DO NOT RUN

$ cat your_file.txt | grep "genomics" >> molecular_grep.txt

$ cat molecular_grep.txt
```

``` bash
# Where are you?
$ pwd

# back to a folder before
$ cd ../ 

# back to the previous directory
$ cd -

# Back several folders
$ cd ../../../ # each ../ is a folder back
```

``` bash
# Downloading reference genome and annotation files
$ mkdir bionfo_course_gdma_2023

$ cd bionfo_course_gdma_2023

$ mkdir reference_genome

$ wget https://ftp.ensembl.org/pub/release-109/fasta/gallus_gallus/dna/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz -O reference_genome/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.dna.toplevel.fa.gz

$ wget https://ftp.ensembl.org/pub/release-109/gtf/gallus_gallus/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.109.gtf.gz -O reference_genome/Gallus_gallus.bGalGal1.mat.broiler.GRCg7b.109.gtf.gz



$ cd reference_genome/
$ gunzip Gallus_gallus*
```
