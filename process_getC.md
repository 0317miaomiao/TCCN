# Computing Random Matrix C(s) Using Graph-Theoretic Approaches

The random matrix **c(s)** can be computed using graph-theoretic approaches. For intuitive illustration, Figure 1 provides a schematic overview of this procedure when exons are of fixed length.

![Computing Process Overview](Image/image_13.png)

## Process Overview

The computation process involves the following five steps:

### Step 1: Identify the Common Regions

Considering a set of transcripts, we divide the transcript space into sections based on the partition of reads and then identify the pre-common regions.

### Step 2: Retain Regions Longer Than Read Length

The pre-common regions identified in Step 1, whose lengths do not exceed the read length, are removed. This ensures that only regions capable of generating reads of the specified length are retained.

### Step 3: Merge Adjacent Exons and Add Edges

After alternative splicing, introns are excised, leaving only exons. Adjacent exons are then merged into a single vertex, with a dotted edge added to connect adjacent vertices. Additionally, a sliding window approach is employed to illustrate the merging process and the addition of dotted edges.

To better demonstrate that this step also holds in the general case (with exons of variable length), we provide an additional example, thereby offering a clearer interpretation of the vertices and edges in this context.

#### Example: Variable Length Exons

Consider a diagram illustrating this step, focusing on a single region as an example. Assume exons e1, e2, e3, and e4 have lengths of **200 bp**, **200 bp**, **300 bp**, and **100 bp**, respectively, while the read length is **250 bp**. A sliding window is applied as presented in Figure 2, starting from the first exon in this region.

![Sliding Window Example](Image/image_14.png)

**Detailed Process:**

1. **First Window Position**: Since exon e1 is shorter than the read length, it cannot independently generate a read. Therefore, e1 and e2 are merged into a single vertex, indicating that reads must span both e1 and e2.

2. **Window Shift**: The sliding window is then moved over a small range, and it is found that the reads can span e1, e2, and e3. A dotted edge is added behind the first vertex to reflect this span.

3. **Second Window Position**: The window is then moved to the second exon, e2. Since e2 is also shorter than the read length, e2 and e3 are merged into a second vertex, indicating that reads must span e2 and e3. However, the sliding window reveals that reads cannot span e2, e3, and e4, so no dotted edge is added behind the second vertex.

4. **Third Window Position**: Next, the window is slid to the third exon, e3, which is longer than the read length and can generate reads on its own. Therefore, e3 and e4 are not merged. However, since the reads can span both e3 and e4, a dotted edge is added behind the third vertex.

5. **Final Window Position**: Finally, the window is moved to the last exon, e4. As e4 is shorter, no further vertices are generated, and no additional dotted edges are added.

### Step 4: Identify Exclusive Regions and Assign Weights

After removing repetitive vertices and edges in all sections, the exclusive regions are identified. This sliding window approach can also be used to interpret the assignment of weights.

**Weight Assignment Example:**
- The **first vertex**, which indicates that reads can span exons e1 and e2, corresponds to **151 possible starting positions** for the reads. Therefore, a weight of **151** is assigned to the first vertex.
- The **dotted edge** behind the first vertex, which indicates that reads can span exons e1, e2, and e3, corresponds to **49 possible starting positions** for the reads, so a weight of **49** is assigned to the dotted edge.
- Similarly, the **second and third vertices** are assigned weights of **200** and **51**, respectively, while the **dotted edge** behind the third vertex is assigned a weight of **100**.

### Step 5: Compute the Random Matrix

The total value in each section is computed by summing the weighted vertices and edges, which represents the **numerator** of the random matrix elements. The **denominator** is given by the length of the respective regions. The elements of the random matrix are then calculated by division.

## Summary

This graph-theoretic approach provides a systematic method for computing the random matrix c(s) by:
1. Identifying relevant genomic regions
2. Modeling transcript structures as graph vertices and edges
3. Assigning appropriate weights based on read mapping possibilities
4. Computing matrix elements through normalization