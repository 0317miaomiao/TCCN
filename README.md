# How to calculate the number of conditions in the transcript space using the relevant code

First through the **process_gene_data** function to process the transcript data, get arrays, gene_dict two data, and then in the **ComputeC** function to input the previously obtained arrays, 
gene_dict data and their own choice of read_length, and finally will get a C matrix. (Specific procedures and code can be through the **Example** of this jupyter notebook to learn.)

For example, for the first example in the code, after processing by the **process_gene_data** function, we can get a list 

[[1, 2, 3, 4, 5, 7, 8, 9, 10, 12, 13, 14, 16, 17, 18, 19, 20, 21, 22, 23, 24], [1, 2, 3, 4, 5, 6], [4, 5, 8, 9, 10, 11], [9, 10, 12, 13, 14, 15], [22, 23, 25, 26]], 

and a dictionary 

{1: 919, 2: 184, 3: 179, 4: 176, 5: 106, 6: 945, 7: 82, 8: 91, 9: 107, 10: 90, 11: 81, 12: 157, 13: 107, 14: 103, 15: 85, 16: 135, 17: 96, 18: 138, 19: 117, 20: 166, 21: 386, 22: 180, 23: 79, 24: 216, 25: 216, 26: 95}. 

This corresponds to the composition of the transcript space of the gene, with the list indicating the position of the exon and the dictionary indicating the length of the corresponding exon. The composition of the transcript space is shown below:

![示例图片](Image/figure_github.png)

After getting the composition of the transcript we use the function ComputeC to calculate the number of conditions under different read lengths, the specific effect is shown in the figure below:

![示例图片](Image/image_12.png)
