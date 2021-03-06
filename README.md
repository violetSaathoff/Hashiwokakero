# Hashiwokakero

## Overview
This project is a demonstration of the fact that the puzzle Hashiwokakero (https://en.wikipedia.org/wiki/Hashiwokakero) is Turing Complete. I demonstrate this by creating circuit elements using Hashiwokakero puzzles. For this to work I need to define an encoding scheme which has inputs/outputs, I need to create bendable wires, I need to be able to split wires, and I need to be able to compute all 1/2 input logical functions (while a NOT gate and an OR gate would suffice, I have puzzles for each standard logic gate: NOT, OR, NOR, AND, NAND, XOR, XNOR, as well as additional puzzles to split wires and cross wires, as well as a non-trivial buffer circuit to protect more complicated circuits).

Data is encoded as edges which have either a single bridge (false) or a double bridge (true). Notably, since (1 + 2) = (2 + 1) = 3, data can be notted using a 3. This also allows for the creationg of bendable wires. For example: 
    
    As a NOT gate (I = input, O = output, F = false, T = true):
    
    I - 3 - O       ||       F - 3 = T       ||       T - 3 - F
    
    As a Bendable Wire:
    
            3 - O       ||               3 - F       ||               3 = T
    I - 3   |           ||       F - 3   ‖           ||       T = 3   |    
        |   |           ||           ‖   ‖           ||           |   |    
        3 - 3           ||           3 - 3           ||           3 = 3    

The left panels are the circuits with forced bridges solved, but arbitrary inputs/outputs. The middle panels are the circuits solved when the input is false. The right panels are the circuits solved when the input is true.

Inputs and outputs are achieved by specifying particular edges as either inputs or outputs. Since each edge has 2 ends, this will always generate 2 copies of each input (and means that each output needs to be calculated twice. In my circuits I position my inputs along the left of the puzzle and my outputs along the right. For example, here is a circuit with 1 input/output that is a unitary logic function (forced bridges are already added):

        2       2           ||           2       2    
        ‖       ‖           ||           ‖       ‖    
    2 = 8 ===== 8 = 2       ||       2 = 8 ===== 8 = 2
        ‖       ‖           ||           ‖       ‖    
    3 - 7 ----- 7 - 3       ||       3 - 7 - 3 - 7 - 3
    |   ‖       ‖   |       ||       |   ‖       ‖   |
    In  6 ===== 6  Out      ||       In  6 ===== 6  Out
    |   ‖       ‖   |       ||       |   ‖       ‖   |
    3 - 7 ----- 7 - 3       ||       3 - 7 - 3 - 7 - 3
        ‖       ‖           ||           ‖       ‖    
    2 = 8 ===== 8 = 2       ||       2 = 8 ===== 8 = 2
        ‖       ‖           ||           ‖       ‖    
        2       2           ||           2       2    
 
The left panel is a trivial buffer, and the right panel is a NOT gate. Note that the 7's act as 3's since the 8's/6's force there to be 4 additional bridges connected to the 7 which don't transmit data, leaving 3 bridges which do connect to the 7 to transmit data. Also, the double-bridge structure is very useful because structures like it can be used to force the puzzle to be fully connected when it otherwise might not be, and it can be used to block communication between vertices which could otherwise talk to each other (an example of blocked communication is the 6's in the puzzle above, which block vertical communication between the 7's). Finally, it's worth noting that the output is effectively a data terminator (the ability to terminate data will be useful later).

The most important building block for making logic gates is what I refer to as a "switch block", or an "if/else block". The basic idea behind this is to leverage the fact that bridges are not allowed to cross to create a structure that controls how data from one input flows based on the truth value of another input. To do this I convert the 2/1 logic levels to 0/1 logic levels using a 2. For example: 

    I - 2   2 - O       ||       F - 2 - 2 - F       ||       T = 2   2 = T

For the time being please ignore the fact that these puzzles are not fully connected (this will be fixed later with the cage structure). Note how the edge between the 2's allows vertical bridges to pass through it when the input is true, but not when the input is false. On it's own, this isn't particularly useful since the bridge which would have passed through the if gate still needs to connect somewhere. To get around this I added a NOT gate followed by a perpendicular if gate (note that since the logic levels are currently 0/1 instead of 2/1, a NOT gate is a 1 instead of a 3).

    I - 2   1               ||       F - 2 - 1               ||       T = 2   1
          X                 ||             X                 ||             X |
            2 - 3 - O       ||               2 = 3 - F       ||               2 - 3 = T
 
In this circuit, the data from X will propagate to the right when the input is false (middle panel), and upwards when the input is true (right panel). Also, note that the input is not destroyed in this process. (Side note: the switch block, the cage, the inputs/outputs, and the terminator can be used to create any arbitrary logical function within my framework! That being said, the Fake XOR makes some tasks much easier.)
 
The easiest circuit where the switch block is useful is a wire splitter. The basic building block of the wire splitter is:
 
          2 --- 3 - O       ||             2 === 3 - F       ||             2 --- 3 = T
    I - 2   1               ||       F - 2 - 1               ||       T = 2 | 1        
          1   2 - O         ||             1 - 2 - F         ||             1 | 2 = T  
            2 - 3 - O       ||               2 = 3 - F       ||               2 - 3 = T
 
Although I often use this as a building block within a larger puzzle, it can be used directly as a wire splitter like:
 
        2           2           ||           2           2           ||           2           2    
        ‖           ‖           ||           ‖           ‖           ||           ‖           ‖    
    2 = 8 ========= 8 = 2       ||       2 = 8 ========= 8 = 2       ||       2 = 8 ========= 8 = 2
        ‖           ‖           ||           ‖           ‖           ||           ‖           ‖    
    3 - 7 --------- 7 - 3       ||       3 = 7 --------- 7 = 3       ||       3 - 7 ========= 7 - 3
    |   ‖ 2 --- 3   ‖   |       ||       |   ‖ 2 --- 3   ‖   |       ||       ‖   ‖ 2 === 3   ‖   ‖
    |   4       |   4   |       ||       |   4 |     ‖   4   |       ||       ‖   4       |   4   ‖
    |   ‖     2 |   ‖   |       ||       |   ‖ |   2 ‖   ‖   |       ||       ‖   ‖     2 |   ‖   ‖
    3 - 6   1 ‖ 3 - 7 - 3       ||       3 = 6 | 1 ‖ 3 - 7 = 3       ||       3 - 6 - 1 ‖ 3 = 7 - 3
        ‖     ‖     ‖           ||           ‖ | | ‖     ‖           ||           ‖     ‖     ‖    
        ‖     4 === 8 = 2       ||           ‖ | | 4 === 8 = 2       ||           ‖     4 === 8 = 2
        ‖           ‖           ||           ‖ | |       ‖           ||           ‖           ‖    
        ‖ 1         6 - 3       ||           ‖ 1 |       6 - 3       ||           ‖ 1 ------- 6 - 3
        ‖           ‖   |       ||           ‖   |       ‖   |       ||           ‖           ‖   ‖
        ‖       2 = 6   |       ||           ‖   |   2 = 6   |       ||           ‖       2 = 6   ‖
        ‖           ‖   |       ||           ‖   |       ‖   |       ||           ‖           ‖   ‖
        ‖   2 ----- 7 - 3       ||           ‖   2 ----- 7 = 3       ||           ‖   2 ===== 7 - 3
        ‖           ‖           ||           ‖           ‖           ||           ‖           ‖    
    2 = 8 ========= 8 = 2       ||       2 = 8 ========= 8 = 2       ||       2 = 8 ========= 8 = 2
        ‖           ‖           ||           ‖           ‖           ||           ‖           ‖    
        2           2           ||           2           2           ||           2           2    

The ability to split wires can also be used to create a non-trivial buffer, which can be used to protect circuits from propagating indeterminate logic:

        2     2           ||           2     2           ||           2     2    
        ‖     ‖           ||           ‖     ‖           ||           ‖     ‖    
    2 = 8 === 8 = 2       ||       2 = 8 === 8 = 2       ||       2 = 8 === 8 = 2
        ‖     ‖           ||           ‖     ‖           ||           ‖     ‖    
        ‖ 2 - 7 - 3       ||           ‖ 2 - 7 = 3       ||           ‖ 2 = 7 - 3
    3 - 6   1 ‖   |       ||       3 = 6 | 1 ‖   |       ||       3 - 6 - 1 ‖   ‖
    |   ‖     4   |       ||       |   ‖ | | 4   |       ||       ‖   ‖     4   ‖
    |   4     ‖   |       ||       |   4 | | ‖   |       ||       ‖   4     ‖   ‖
    |   ‖ 1   6 - 3       ||       |   ‖ 1 | 6 = 3       ||       ‖   ‖ 1 - 6 - 3
    3 - 7 - 2 ‖           ||       3 = 7 - 2 ‖           ||       3 - 7 = 2 ‖    
        ‖     ‖           ||           ‖     ‖           ||           ‖     ‖    
    2 = 8 === 8 = 2       ||       2 = 8 === 8 = 2       ||       2 = 8 === 8 = 2
        ‖     ‖           ||           ‖     ‖           ||           ‖     ‖    
        2     2           ||           2     2           ||           2     2    

Since the input is terminated, the 2 halves of the input are foced be have the same truth value. Furthermore, since the output is generated using a loop, both halves of the output are also forced to have the same truth value. This isn't important here, but in more complicated circuits it's potentially possible for there to be circuits which allow extra solutions under certain input combinations that lead to ambiguous logic, but buffers can prevent this from being an issue.

Although this isn't actually necessary for Turing Competeness, the next circuit element I want to explain is the FAKE XOR. The Fake XOR is works by having a vertex that has 2 inputs (A, B) and 1 output (C). If the value on the vertex is V, then we can compute C using the equation: C = V - A - B. Notably, this equation preserves parity (which is why this gate is most akin to an xor gate). Unfortunately, for the resulting data to be a valid logic level we need there to be only 2 possible values of C. This means that at best there are 3 possible input combinations which satisfy the gate (this is why it's a FAKE XOR and not a real XOR). For example, if V = 5, then A and B can't simultaneously be 1 since then C would have to be 3, but that's illegal (other options work: 5 - 2 - 2 = 1, 5 - 2 - 1 = 2, 5 - 1 - 2 = 2). Alternatively, if V = 4, then A and B can't simultaneously be 2 since then C would have to be 0, which isn't a valid logic level in my 1/2 encoding scheme (other options work: 4 - 1 - 1 = 2, 4 - 1 - 2 = 1, 4 - 2 - 1 = 1). 

Here is a FAKE XOR which uses a 4 (each panel is for a different input combination, and is fully solved when possible - note that because NOT gates are cheap I don't really care that this is actually a FAKE XNOR instead of a FAKE XOR):

    F           ||       F           ||       T           ||       T    
    |           ||       |           ||       ‖           ||       ‖    
    4 = T       ||       4 - F       ||       4 - F       ||       4   X
    |           ||       ‖           ||       |           ||       ‖    
    F           ||       T           ||       F           ||       T    

And here is a FAKE XOR which uses a 5:

    F           ||       F           ||       T           ||       T    
    |           ||       |           ||       ‖           ||       ‖    
    5 = X       ||       5 = T       ||       5 = T       ||       5 - F
    |           ||       ‖           ||       |           ||       ‖    
    F           ||       T           ||       F           ||       T    

Using an IF/ELSE block and 2 FAKE XORs it is actually possible to make an XOR gate which operates on individual data channels!

        2                 2           ||           2                 2           ||           2                 2           ||           2                 2    
        ‖                 ‖           ||           ‖                 ‖           ||           ‖                 ‖           ||           ‖                 ‖    
    2 = 8 =============== 8 = 2       ||       2 = 8 =============== 8 = 2       ||       2 = 8 =============== 8 = 2       ||       2 = 8 =============== 8 = 2
        ‖   2 ======= 3   ‖           ||           ‖   2 ======= 3   ‖           ||           ‖   2 ------- 3   ‖           ||           ‖   2 ======= 3   ‖    
    F - 6 --- 1       |   ‖           ||       F - 6 --- 1       |   ‖           ||       T = 6   | 1       ‖   ‖           ||       T = 6     1       |   ‖    
        ‖   2 --- 2   |   ‖           ||           ‖   2     2   |   ‖           ||           ‖   2 |   2   ‖   ‖           ||           ‖   2 |   2   |   ‖    
    2 = 6   |     |   4 = 7 - F       ||       2 = 6   ‖     ‖   4 - 7 = T       ||       2 = 6   | |   ‖   4 - 7 = T       ||       2 = 6   ‖ |   ‖   4 = 7 - F
        ‖   |     |   |   ‖           ||           ‖   ‖     ‖   ‖   ‖           ||           ‖   | |   ‖   |   ‖           ||           ‖   ‖ |   ‖   |   ‖    
    F - 7 = 3     |   |   ‖           ||       T = 7 - 3     ‖   ‖   ‖           ||       F - 7 = 3 |   ‖   |   ‖           ||       T = 7 - 3 |   ‖   |   ‖    
        ‖     2 = 5 = 3   ‖           ||           ‖     2 = 5 - 3   ‖           ||           ‖     2 - 5 = 3   ‖           ||           ‖     2 - 5 = 3   ‖    
    2 = 8 =============== 8 = 2       ||       2 = 8 =============== 8 = 2       ||       2 = 8 =============== 8 = 2       ||       2 = 8 =============== 8 = 2
        ‖                 ‖           ||           ‖                 ‖           ||           ‖                 ‖           ||           ‖                 ‖    
        2                 2           ||           2                 2           ||           2                 2           ||           2                 2    

Next, using this XOR element and the SWAP element it is possible to create an actual XOR gate (note that the output is buffered)! There's actually a trick here where I merged the SWAP element with one of the XOR elements, which is why one of the FAKE XORs actually uses a 3 (at position (4,4)).

        2                     2           ||           2                     2           ||           2                     2           ||           2                     2    
        ‖                     ‖           ||           ‖                     ‖           ||           ‖                     ‖           ||           ‖                     ‖    
    2 = 8 =================== 8 = 2       ||       2 = 8 =================== 8 = 2       ||       2 = 8 =================== 8 = 2       ||       2 = 8 =================== 8 = 2
        ‖                     ‖           ||           ‖                     ‖           ||           ‖                     ‖           ||           ‖                     ‖    
    3 = 7 - 3 - 5 = 4 = 4 --- 7 = 3       ||       3 = 7 - 3 = 5 = 4 = 4 === 7 - 3       ||       3 - 7 = 3 - 5 = 4 = 4 === 7 - 3       ||       3 - 7 = 3 - 5 = 4 = 4 --- 7 = 3
    |   ‖   | 2 ‖       |     ‖   |       ||       |   ‖     2 |             ‖   ‖       ||       ‖   ‖     2 ‖             ‖   ‖       ||       ‖   ‖     2 ‖       |     ‖   |
    |   4   | ‖ 4       | 3 = 6   |       ||       |   4     ‖ 4 ------- 3 = 6   ‖       ||       ‖   4     ‖ 4 ------- 3 = 6   ‖       ||       ‖   4     ‖ 4       | 3 = 6   |
    |   ‖   | ‖ ‖       | |   ‖   |       ||       |   ‖     ‖ ‖             ‖   ‖       ||       ‖   ‖     ‖ |             ‖   ‖       ||       ‖   ‖     ‖ ‖       | |   ‖   |
    3 = 6   | 3 ‖       1 |   6 = 3       ||       3 = 6     3 ‖       1 --- 6 - 3       ||       3 - 6 --- 3 |       1 --- 6 - 3       ||       3 - 6 --- 3 ‖       1 |   6 = 3
        ‖   | | ‖         |   ‖           ||           ‖     | ‖             ‖           ||           ‖       |             ‖           ||           ‖       ‖         |   ‖    
    2 = 6   | | ‖ 2 = 5 - 6 = 8 = 2       ||       2 = 6     | ‖ 2 - 5 = 6 = 8 = 2       ||       2 = 6       | 2 = 5 = 6 = 8 = 2       ||       2 = 6       ‖ 2 - 5 = 6 = 8 = 2
        ‖   | | ‖     ‖   ‖   ‖           ||           ‖     | ‖ |   ‖   ‖   ‖           ||           ‖       |     |   ‖   ‖           ||           ‖       ‖ |   ‖   |   ‖    
    3 = 7 - 2 | 2     ‖   3   ‖           ||       3 - 7 = 2 | 2 |   ‖   3   ‖           ||       3 = 7 - 2 - 2     |   3   ‖           ||       3 - 7 = 2   2 |   ‖   3   ‖    
    |   ‖     |       ‖   |   ‖           ||       ‖   ‖     |   |   ‖   |   ‖           ||       |   ‖             |   |   ‖           ||       ‖   ‖         |   ‖   ‖   ‖    
    |   6 === 3       2   |   ‖           ||       ‖   6 === 3   |   2   |   ‖           ||       |   6 === 3 ----- 2   |   ‖           ||       ‖   6 === 3   |   2   ‖   ‖    
    |   ‖                 |   ‖           ||       ‖   ‖         |       |   ‖           ||       |   ‖                 |   ‖           ||       ‖   ‖     |   |       ‖   ‖    
    3 = 7 - 2 --- 1       |   ‖           ||       3 - 7 = 2     1       |   ‖           ||       3 = 7 - 2 --- 1       |   ‖           ||       3 - 7 = 2 |   1       ‖   ‖    
        ‖     2 ========= 3   ‖           ||           ‖     2 ========= 3   ‖           ||           ‖     2 ========= 3   ‖           ||           ‖     2 --------- 3   ‖    
    2 = 8 =================== 8 = 2       ||       2 = 8 =================== 8 = 2       ||       2 = 8 =================== 8 = 2       ||       2 = 8 =================== 8 = 2
        ‖                     ‖           ||           ‖                     ‖           ||           ‖                     ‖           ||           ‖                     ‖    
        2                     2           ||           2                     2           ||           2                     2           ||           2                     2    


Finally, to be Turing Complete it has to be possible to compute an OR logic Function (all panels are fully solved): 

        2                         2           ||           2                         2           ||           2                         2           ||           2                         2    
        ‖                         ‖           ||           ‖                         ‖           ||           ‖                         ‖           ||           ‖                         ‖    
    2 = 8 ======================= 8 = 2       ||       2 = 8 ======================= 8 = 2       ||       2 = 8 ======================= 8 = 2       ||       2 = 8 ======================= 8 = 2
        ‖                         ‖           ||           ‖                         ‖           ||           ‖                         ‖           ||           ‖                         ‖    
    3 = 7 ------------- 3         ‖           ||       3 = 7 ------------- 3         ‖           ||       3 - 7 ============= 3         ‖           ||       3 - 7 ============= 3         ‖    
    |   ‖               ‖         ‖           ||       |   ‖               ‖         ‖           ||       ‖   ‖               |         ‖           ||       ‖   ‖               |         ‖    
    |   6 === 4 = 2     ‖         ‖           ||       |   6 === 4 = 2     ‖         ‖           ||       ‖   6 === 4 - 2     |         ‖           ||       ‖   6 === 4 = 2     |         ‖    
    |   ‖               ‖         ‖           ||       |   ‖               ‖         ‖           ||       ‖   ‖     |   |     |         ‖           ||       ‖   ‖               |         ‖    
    3 = 7 - 2 ----- 1   ‖         ‖           ||       3 = 7 - 2 ----- 1   ‖         ‖           ||       3 - 7 = 2 |   | 1   |         ‖           ||       3 - 7 = 2       1   |         ‖    
        ‖               2         ‖           ||           ‖               2         ‖           ||           ‖     |   | |   2         ‖           ||           ‖           |   2         ‖    
    2 = 8 === 3 --------- 2 - 5 = 8 = 2       ||       2 = 8 === 3           2 = 5 = 8 = 2       ||       2 = 8 === 3   | |   | 2 = 5 = 8 = 2       ||       2 = 8 === 3     |   | 2 = 5 = 8 = 2
        ‖         2 ------- 2 ‖   ‖           ||           ‖     |   2         2 |   ‖           ||           ‖         2 |   |   2 |   ‖           ||           ‖     |   2 |   |   2 |   ‖    
    3 = 7 - 2 - 1 | 2 = 4   | 3 - 7 = 3       ||       3 - 7 = 2 | 1 ‖ 2 = 4   ‖ 3 = 7 - 3       ||       3 = 7 - 2 - 1 | 2 - 4   ‖ 3 = 7 - 3       ||       3 - 7 = 2 | 1 ‖ 2 - 4   ‖ 3 = 7 - 3
    |   ‖         |     ‖   |     ‖   |       ||       ‖   ‖     | | ‖     ‖   ‖     ‖   ‖       ||       |   ‖         |     ‖   ‖     ‖   ‖       ||       ‖   ‖     | | ‖     ‖   ‖     ‖   ‖
    |   6 === 3 - 4 === 8 = 7 === 6   |       ||       ‖   6 === 3 | 4 === 8 = 7 === 6   ‖       ||       |   6 === 3 - 4 === 8 = 7 === 6   ‖       ||       ‖   6 === 3 | 4 === 8 = 7 === 6   ‖
    |   ‖               ‖   ‖     ‖   |       ||       ‖   ‖       |       ‖   |     ‖   ‖       ||       |   ‖               ‖   |     ‖   ‖       ||       ‖   ‖       |       ‖   |     ‖   ‖
    3 = 7 - 3 = 2       2   3 --- 7 = 3       ||       3 - 7 = 3 - 2       2   3 === 7 - 3       ||       3 = 7 - 3 = 2       2   3 === 7 - 3       ||       3 - 7 = 3 - 2       2   3 === 7 - 3
        ‖                         ‖           ||           ‖                         ‖           ||           ‖                         ‖           ||           ‖                         ‖    
    2 = 8 ======================= 8 = 2       ||       2 = 8 ======================= 8 = 2       ||       2 = 8 ======================= 8 = 2       ||       2 = 8 ======================= 8 = 2
        ‖                         ‖           ||           ‖                         ‖           ||           ‖                         ‖           ||           ‖                         ‖    
        2                         2           ||           2                         2           ||           2                         2           ||           2                         2    

The OR gate works by using a doubled IF/ELSE block controlled by the first input (A) to only allow the second input (B) to talk to the output when the A is false. Meanwhile, if the B is false it will send bridges through whichever gate is open (the upper/IF gate when A is true which leads to a terminator, and the right/ELSE gate when A is false which leads to the output). This means that the only time bridges pass through the ELSE gate to change the output is when both inputs are false. When this happens the output is false, otherwise the output is true, which is the desired logic function. Also, note that both inputs are buffered (A is buffered because its input channels happend to loop together, and B is buffered using an actual buffer circuit).

The NOR, AND, and NAND gates are all made by notting the inputs/outputs of the OR gate appropriately; for completeness I've included their solved states below. First is the NOR gate:

        2                       2           ||           2                       2           ||           2                       2           ||           2                       2    
        ‖                       ‖           ||           ‖                       ‖           ||           ‖                       ‖           ||           ‖                       ‖    
    2 = 8 ===================== 8 = 2       ||       2 = 8 ===================== 8 = 2       ||       2 = 8 ===================== 8 = 2       ||       2 = 8 ===================== 8 = 2
        ‖                       ‖           ||           ‖                       ‖           ||           ‖                       ‖           ||           ‖                       ‖    
    3 = 7 ------------- 3       ‖           ||       3 = 7 ------------- 3       ‖           ||       3 - 7 ============= 3       ‖           ||       3 - 7 ============= 3       ‖    
    |   ‖               ‖       ‖           ||       |   ‖               ‖       ‖           ||       ‖   ‖               |       ‖           ||       ‖   ‖               |       ‖    
    |   6 === 4 = 2     ‖       ‖           ||       |   6 === 4 = 2     ‖       ‖           ||       ‖   6 === 4 - 2     |       ‖           ||       ‖   6 === 4 = 2     |       ‖    
    |   ‖               ‖       ‖           ||       |   ‖               ‖       ‖           ||       ‖   ‖     |   |     |       ‖           ||       ‖   ‖               |       ‖    
    3 = 7 - 2 ----- 1   ‖       ‖           ||       3 = 7 - 2 ----- 1   ‖       ‖           ||       3 - 7 = 2 |   | 1   |       ‖           ||       3 - 7 = 2       1   |       ‖    
        ‖               2       ‖           ||           ‖               2       ‖           ||           ‖     |   | |   2       ‖           ||           ‖           |   2       ‖    
    2 = 8 === 3 ----------- 4 = 8 = 2       ||       2 = 8 === 3             4 = 8 = 2       ||       2 = 8 === 3   | |   |   4 = 8 = 2       ||       2 = 8 === 3     |   |   4 = 8 = 2
        ‖         2 ----- 2 |   ‖           ||           ‖     |   2       2 ‖   ‖           ||           ‖         2 |   | 2 ‖   ‖           ||           ‖     |   2 |   | 2 ‖   ‖    
    3 = 7 - 2 - 1 | 2 = 4 | 3 = 7 - 3       ||       3 - 7 = 2 | 1 ‖ 2 = 4 ‖ 3 - 7 = 3       ||       3 = 7 - 2 - 1 | 2 - 4 ‖ 3 - 7 = 3       ||       3 - 7 = 2 | 1 ‖ 2 - 4 ‖ 3 - 7 = 3
    |   ‖         |     ‖ |     ‖   ‖       ||       ‖   ‖     | | ‖     ‖ ‖     ‖   |       ||       |   ‖         |     ‖ ‖     ‖   |       ||       ‖   ‖     | | ‖     ‖ ‖     ‖   |
    |   6 === 3 - 4 === 6 | 2 = 6   ‖       ||       ‖   6 === 3 | 4 === 6 ‖ 2 = 6   |       ||       |   6 === 3 - 4 === 6 ‖ 2 = 6   |       ||       ‖   6 === 3 | 4 === 6 ‖ 2 = 6   |
    |   ‖               ‖ |     ‖   ‖       ||       ‖   ‖       |       ‖ ‖     ‖   |       ||       |   ‖               ‖ ‖     ‖   |       ||       ‖   ‖       |       ‖ ‖     ‖   |
    3 = 7 - 3 = 2       ‖ 3 === 7 - 3       ||       3 - 7 = 3 - 2       ‖ 3 --- 7 = 3       ||       3 = 7 - 3 = 2       ‖ 3 --- 7 = 3       ||       3 - 7 = 3 - 2       ‖ 3 --- 7 = 3
        ‖               ‖       ‖           ||           ‖               ‖       ‖           ||           ‖               ‖       ‖           ||           ‖               ‖       ‖    
    2 = 8 ============= 6 ===== 8 = 2       ||       2 = 8 ============= 6 ===== 8 = 2       ||       2 = 8 ============= 6 ===== 8 = 2       ||       2 = 8 ============= 6 ===== 8 = 2
        ‖                       ‖           ||           ‖                       ‖           ||           ‖                       ‖           ||           ‖                       ‖    
        2                       2           ||           2                       2           ||           2                       2           ||           2                       2    

Next is the AND gate:

        2                       2           ||           2                       2           ||           2                       2           ||           2                       2    
        ‖                       ‖           ||           ‖                       ‖           ||           ‖                       ‖           ||           ‖                       ‖    
    2 = 8 ===================== 8 = 2       ||       2 = 8 ===================== 8 = 2       ||       2 = 8 ===================== 8 = 2       ||       2 = 8 ===================== 8 = 2
        ‖                       ‖           ||           ‖                       ‖           ||           ‖                       ‖           ||           ‖                       ‖    
    3 = 7 ------------- 2       ‖           ||       3 = 7 ------------- 2       ‖           ||       3 - 7 ============= 2       ‖           ||       3 - 7 ============= 2       ‖    
    |   ‖               |       ‖           ||       |   ‖               |       ‖           ||       ‖   ‖                       ‖           ||       ‖   ‖                       ‖    
    |   6 === 4 = 2     |       ‖           ||       |   6 === 4 - 2     |       ‖           ||       ‖   6 === 4 = 2             ‖           ||       ‖   6 === 4 = 2             ‖    
    |   ‖               |       ‖           ||       |   ‖     |   |     |       ‖           ||       ‖   ‖                       ‖           ||       ‖   ‖                       ‖    
    3 = 6           1   |       ‖           ||       3 = 6     |   | 1   |       ‖           ||       3 - 6 --------- 1           ‖           ||       3 - 6 --------- 1           ‖    
        ‖           |   |       ‖           ||           ‖     |   | |   |       ‖           ||           ‖                       ‖           ||           ‖                       ‖    
    2 = 8 === 3     |   |   4 = 8 = 2       ||       2 = 8 === 3   | |   |   4 = 8 = 2       ||       2 = 8 === 3             4 = 8 = 2       ||       2 = 8 === 3 ----------- 4 = 8 = 2
        ‖     |   2 |   | 2 ‖   ‖           ||           ‖         2 |   | 2 ‖   ‖           ||           ‖     |   2       2 ‖   ‖           ||           ‖         2 ----- 2 |   ‖    
    3 = 6     | 1 ‖ 2 - 4 ‖ 3 - 7 = 3       ||       3 - 6 ----- 1 | 2 - 4 ‖ 3 - 7 = 3       ||       3 = 6     | 1 ‖ 2 = 4 ‖ 3 - 7 = 3       ||       3 - 6 ----- 1 | 2 = 4 | 3 = 7 - 3
    |   ‖     | | ‖     ‖ ‖     ‖   |       ||       ‖   ‖         |     ‖ ‖     ‖   |       ||       |   ‖     | | ‖     ‖ ‖     ‖   |       ||       ‖   ‖         |     ‖ |     ‖   ‖
    |   6 === 3 | 4 === 6 ‖ 2 = 6   |       ||       ‖   6 === 3 - 4 === 6 ‖ 2 = 6   |       ||       |   6 === 3 | 4 === 6 ‖ 2 = 6   |       ||       ‖   6 === 3 - 4 === 6 | 2 = 6   ‖
    |   ‖       |       ‖ ‖     ‖   |       ||       ‖   ‖               ‖ ‖     ‖   |       ||       |   ‖       |       ‖ ‖     ‖   |       ||       ‖   ‖               ‖ |     ‖   ‖
    3 = 7 ----- 2       ‖ 3 --- 7 = 3       ||       3 - 7 ===== 2       ‖ 3 --- 7 = 3       ||       3 = 7 ----- 2       ‖ 3 --- 7 = 3       ||       3 - 7 ===== 2       ‖ 3 === 7 - 3
        ‖               ‖       ‖           ||           ‖               ‖       ‖           ||           ‖               ‖       ‖           ||           ‖               ‖       ‖    
    2 = 8 ============= 6 ===== 8 = 2       ||       2 = 8 ============= 6 ===== 8 = 2       ||       2 = 8 ============= 6 ===== 8 = 2       ||       2 = 8 ============= 6 ===== 8 = 2
        ‖                       ‖           ||           ‖                       ‖           ||           ‖                       ‖           ||           ‖                       ‖    
        2                       2           ||           2                       2           ||           2                       2           ||           2                       2    

Next, the NAND gate:

        2                         2           ||           2                         2           ||           2                         2           ||           2                         2    
        ‖                         ‖           ||           ‖                         ‖           ||           ‖                         ‖           ||           ‖                         ‖    
    2 = 8 ======================= 8 = 2       ||       2 = 8 ======================= 8 = 2       ||       2 = 8 ======================= 8 = 2       ||       2 = 8 ======================= 8 = 2
        ‖                         ‖           ||           ‖                         ‖           ||           ‖                         ‖           ||           ‖                         ‖    
    3 = 7 ------------- 2         ‖           ||       3 = 7 ------------- 2         ‖           ||       3 - 7 ============= 2         ‖           ||       3 - 7 ============= 2         ‖    
    |   ‖               |         ‖           ||       |   ‖               |         ‖           ||       ‖   ‖                         ‖           ||       ‖   ‖                         ‖    
    |   6 === 4 = 2     |         ‖           ||       |   6 === 4 - 2     |         ‖           ||       ‖   6 === 4 = 2               ‖           ||       ‖   6 === 4 = 2               ‖    
    |   ‖               |         ‖           ||       |   ‖     |   |     |         ‖           ||       ‖   ‖                         ‖           ||       ‖   ‖                         ‖    
    3 = 6           1   |         ‖           ||       3 = 6     |   | 1   |         ‖           ||       3 - 6 --------- 1             ‖           ||       3 - 6 --------- 1             ‖    
        ‖           |   |         ‖           ||           ‖     |   | |   |         ‖           ||           ‖                         ‖           ||           ‖                         ‖    
    2 = 8 === 3     |   | 2 = 5 = 8 = 2       ||       2 = 8 === 3   | |   | 2 = 5 = 8 = 2       ||       2 = 8 === 3           2 = 5 = 8 = 2       ||       2 = 8 === 3 --------- 2 - 5 = 8 = 2
        ‖     |   2 |   |   2 |   ‖           ||           ‖         2 |   |   2 |   ‖           ||           ‖     |   2         2 |   ‖           ||           ‖         2 ------- 2 ‖   ‖    
    3 = 6     | 1 ‖ 2 - 4   ‖ 3 = 7 - 3       ||       3 - 6 ----- 1 | 2 - 4   ‖ 3 = 7 - 3       ||       3 = 6     | 1 ‖ 2 = 4   ‖ 3 = 7 - 3       ||       3 - 6 ----- 1 | 2 = 4   | 3 - 7 = 3
    |   ‖     | | ‖     ‖   ‖     ‖   ‖       ||       ‖   ‖         |     ‖   ‖     ‖   ‖       ||       |   ‖     | | ‖     ‖   ‖     ‖   ‖       ||       ‖   ‖         |     ‖   |     ‖   |
    |   6 === 3 | 4 === 8 = 7 === 6   ‖       ||       ‖   6 === 3 - 4 === 8 = 7 === 6   ‖       ||       |   6 === 3 | 4 === 8 = 7 === 6   ‖       ||       ‖   6 === 3 - 4 === 8 = 7 === 6   |
    |   ‖       |       ‖   |     ‖   ‖       ||       ‖   ‖               ‖   |     ‖   ‖       ||       |   ‖       |       ‖   |     ‖   ‖       ||       ‖   ‖               ‖   ‖     ‖   |
    3 = 7 ----- 2       2   3 === 7 - 3       ||       3 - 7 ===== 2       2   3 === 7 - 3       ||       3 = 7 ----- 2       2   3 === 7 - 3       ||       3 - 7 ===== 2       2   3 --- 7 = 3
        ‖                         ‖           ||           ‖                         ‖           ||           ‖                         ‖           ||           ‖                         ‖    
    2 = 8 ======================= 8 = 2       ||       2 = 8 ======================= 8 = 2       ||       2 = 8 ======================= 8 = 2       ||       2 = 8 ======================= 8 = 2
        ‖                         ‖           ||           ‖                         ‖           ||           ‖                         ‖           ||           ‖                         ‖    
        2                         2           ||           2                         2           ||           2                         2           ||           2                         2    

And lastly, the XNOR gate, which is an XOR gate with a couple NOT gates deleted:

        2                     2           ||           2                     2           ||           2                     2           ||           2                     2    
        ‖                     ‖           ||           ‖                     ‖           ||           ‖                     ‖           ||           ‖                     ‖    
    2 = 8 =================== 8 = 2       ||       2 = 8 =================== 8 = 2       ||       2 = 8 =================== 8 = 2       ||       2 = 8 =================== 8 = 2
        ‖                     ‖           ||           ‖                     ‖           ||           ‖                     ‖           ||           ‖                     ‖    
    3 = 7 - 3 - 5 = 6 = 4 === 7 - 3       ||       3 = 7 - 3 = 5 = 6 = 4 --- 7 = 3       ||       3 - 7 = 3 - 5 = 6 = 4 --- 7 = 3       ||       3 - 7 = 3 - 5 = 6 = 4 === 7 - 3
    |   ‖   | 2 ‖   ‖         ‖   ‖       ||       |   ‖     2 |   ‖   |     ‖   |       ||       ‖   ‖     2 ‖   ‖   |     ‖   |       ||       ‖   ‖     2 ‖   ‖         ‖   ‖
    |   4   | ‖ 5 - 4 --- 3 = 6   ‖       ||       |   4     ‖ 5 = 4   | 3 = 6   |       ||       ‖   4     ‖ 5 = 4   | 3 = 6   |       ||       ‖   4     ‖ 5 - 4 --- 3 = 6   ‖
    |   ‖   | ‖ ‖             ‖   ‖       ||       |   ‖     ‖ ‖       | |   ‖   |       ||       ‖   ‖     ‖ |       | |   ‖   |       ||       ‖   ‖     ‖ ‖             ‖   ‖
    3 = 6   | 3 ‖       1 --- 6 - 3       ||       3 = 6     3 ‖       1 |   6 = 3       ||       3 - 6 --- 3 |       1 |   6 = 3       ||       3 - 6 --- 3 ‖       1 --- 6 - 3
        ‖   | | ‖             ‖           ||           ‖     | ‖         |   ‖           ||           ‖       |         |   ‖           ||           ‖       ‖             ‖    
    2 = 6   | | ‖ 2 - 5 = 6 = 8 = 2       ||       2 = 6     | ‖ 2 = 5 - 6 = 8 = 2       ||       2 = 6       | 2 - 5 = 6 = 8 = 2       ||       2 = 6       ‖ 2 = 5 = 6 = 8 = 2
        ‖   | | ‖ |   ‖   ‖   ‖           ||           ‖     | ‖     ‖   ‖   ‖           ||           ‖       | |   ‖   |   ‖           ||           ‖       ‖     |   ‖   ‖    
    3 = 7 - 2 | 2 |   ‖   3   ‖           ||       3 - 7 = 2 | 2     ‖   3   ‖           ||       3 = 7 - 2 - 2 |   ‖   3   ‖           ||       3 - 7 = 2   2     |   3   ‖    
    |   ‖     |   |   ‖   |   ‖           ||       ‖   ‖     |       ‖   |   ‖           ||       |   ‖         |   ‖   ‖   ‖           ||       ‖   ‖             |   |   ‖    
    |   6 === 3   |   2   |   ‖           ||       ‖   6 === 3       2   |   ‖           ||       |   6 === 3   |   2   ‖   ‖           ||       ‖   6 === 3 ----- 2   |   ‖    
    |   ‖         |       |   ‖           ||       ‖   ‖                 |   ‖           ||       |   ‖     |   |       ‖   ‖           ||       ‖   ‖                 |   ‖    
    3 = 6         1       |   ‖           ||       3 - 6 ------- 1       |   ‖           ||       3 = 6     |   1       ‖   ‖           ||       3 - 6 ------- 1       |   ‖    
        ‖     2 ========= 3   ‖           ||           ‖     2 ========= 3   ‖           ||           ‖     2 --------- 3   ‖           ||           ‖     2 ========= 3   ‖    
    2 = 8 =================== 8 = 2       ||       2 = 8 =================== 8 = 2       ||       2 = 8 =================== 8 = 2       ||       2 = 8 =================== 8 = 2
        ‖                     ‖           ||           ‖                     ‖           ||           ‖                     ‖           ||           ‖                     ‖    
        2                     2           ||           2                     2           ||           2                     2           ||           2                     2    

Finally, although it is possible to swap elements using an xor swap, in practice there is are better options which are made of an IF/ELSE gate and a FAKE XOR gate. The first option is:

          2 = 5 --- F       ||             2 = 5 === T       ||             2 = 5 --- F       ||             2 - 5 === T
    F - 2 - 1 ‖             ||       F - 2 - 1 |             ||       T = 2   1 ‖             ||       T = 2 | 1 ‖      
          2   2             ||             2 - 2             ||             2 | 2             ||             2 | 2      
          ‖ 2 = 3 - F       ||             | 2 = 3 - F       ||             ‖ 2 - 3 = T       ||             | 2 - 3 = T
    F --- 3                 ||       T === 3                 ||       F --- 5 = 2             ||       T === 5 = 2      

The second option is:

          2 = 3             ||             2 = 3             ||             2 - 3             ||             2 = 3      
    F - 2 - 1 |             ||       F - 2 - 1 |             ||       T = 2 | 1 ‖             ||       T = 2   1 |      
    F --- 2 - 3 --- F       ||       T === 2   3 === T       ||       F --- 2 | 3 --- F       ||       T === 2 | 3 === T
            2 = 3 - F       ||               2 = 3 - F       ||               2 - 3 = T       ||               2 - 3 = T

Next, although there are several ways to use this element to make an actual SWAP gate, the best option I've found is to pair the IF/ELSE gates of the first input while splitting the second input:

        2                               2           ||           2                               2           ||           2                               2           ||           2                               2    
        ‖                               ‖           ||           ‖                               ‖           ||           ‖                               ‖           ||           ‖                               ‖    
    2 = 8 === 6 ======================= 8 = 2       ||       2 = 8 === 6 ======================= 8 = 2       ||       2 = 8 === 6 ======================= 8 = 2       ||       2 = 8 === 6 ======================= 8 = 2
        ‖     ‖                         ‖           ||           ‖     ‖                         ‖           ||           ‖     ‖                         ‖           ||           ‖     ‖                         ‖    
    3 = 7 - 3 ‖     3 ============= 3 - 7 = 3       ||       3 = 7 - 3 ‖     3 ------------- 3 = 7 - 3       ||       3 - 7 = 3 ‖     3 ============= 3 - 7 = 3       ||       3 - 7 = 3 ‖     3 ------------- 3 = 7 - 3
    |   ‖   ‖ ‖     | 4 ===== 4         ‖   |       ||       |   ‖   ‖ ‖     ‖ 4 ===== 4         ‖   ‖       ||       ‖   ‖   | ‖     | 4 ===== 4         ‖   |       ||       ‖   ‖   | ‖     ‖ 4 ===== 4         ‖   ‖
    |   4   ‖ 4 === 5 ‖       ‖     2 = 6   |       ||       |   4   ‖ 4 --- 5 ‖       ‖     2 = 6   ‖       ||       ‖   4   | 4 === 5 ‖       ‖     2 = 6   |       ||       ‖   4   | 4 === 5 ‖       ‖     2 = 6   ‖
    |   ‖   2     1 ‖ ‖       ‖         ‖   |       ||       |   ‖   2 |   1 ‖ ‖       ‖         ‖   ‖       ||       ‖   ‖   2 --- 1 ‖ ‖       ‖         ‖   |       ||       ‖   ‖   2 --- 1 | ‖       ‖         ‖   ‖
    3 = 6       1 | ‖ ‖ 2 === 7 === 3 - 7 = 3       ||       3 = 6     | 1 | ‖ ‖ 2 === 7 --- 3 = 7 - 3       ||       3 - 6 ----- 1   ‖ ‖ 2 --- 7 === 3 - 7 = 3       ||       3 - 6 ----- 1   | ‖ 2 === 7 --- 3 = 7 - 3
        ‖     2 | | 2 ‖       |         ‖           ||           ‖     2 | | 2 ‖       ‖         ‖           ||           ‖     2     2 ‖ |     ‖         ‖           ||           ‖     2 --- 2 ‖       ‖         ‖    
    2 = 8 = 2 ‖ | 2 - 4 --- 1 | 2 ===== 8 = 2       ||       2 = 8 = 2 | | 2 - 4 --- 1 ‖ 2 ===== 8 = 2       ||       2 = 8 = 2 ‖   2 = 4 |   1 ‖ 2 ===== 8 = 2       ||       2 = 8 = 2 |   2 = 4     1 ‖ 2 ===== 8 = 2
        ‖     ‖ 4 - 2 --- 1   |         ‖           ||           ‖     | 4 - 2 --- 1   ‖         ‖           ||           ‖     ‖ 4 = 2   | 1 | ‖         ‖           ||           ‖     | 4 = 2     1 | ‖         ‖    
    3 = 7 --- 3 ‖             | 3 = 3 - 7 = 3       ||       3 - 7 === 3 ‖             ‖ 3 = 3 - 7 = 3       ||       3 = 7 --- 3 ‖       | | | ‖ 3 - 3 = 7 - 3       ||       3 - 7 === 3 ‖         | | ‖ 3 - 3 = 7 - 3
    |   ‖       ‖ 3 --- 4 --- 2 |       ‖   |       ||       ‖   ‖       ‖ 3 === 4     2 |       ‖   |       ||       |   ‖       ‖ 3 --- 4 | | 2 ‖       ‖   ‖       ||       ‖   ‖       ‖ 3 === 4 | | 2 ‖       ‖   ‖
    |   6 ===== 4 ‖     ‖   2 = 5 ===== 6   |       ||       ‖   6 ===== 4 |     ‖   2 = 5 ===== 6   |       ||       |   6 ===== 4 ‖     ‖ | 2 - 5 ===== 6   ‖       ||       ‖   6 ===== 4 |     ‖ | 2 - 5 ===== 6   ‖
    |   ‖         ‖     ‖               ‖   |       ||       ‖   ‖         |     ‖               ‖   |       ||       |   ‖         ‖     ‖ |             ‖   ‖       ||       ‖   ‖         |     ‖ |             ‖   ‖
    3 = 7 ------- 3     ‖ 2 ======= 3 - 7 = 3       ||       3 - 7 ======= 3     ‖ 2 ======= 3 - 7 = 3       ||       3 = 7 ------- 3     ‖ 2 ------- 3 = 7 - 3       ||       3 - 7 ======= 3     ‖ 2 ------- 3 = 7 - 3
        ‖               ‖               ‖           ||           ‖               ‖               ‖           ||           ‖               ‖               ‖           ||           ‖               ‖               ‖    
    2 = 8 ============= 6 ============= 8 = 2       ||       2 = 8 ============= 6 ============= 8 = 2       ||       2 = 8 ============= 6 ============= 8 = 2       ||       2 = 8 ============= 6 ============= 8 = 2
        ‖                               ‖           ||           ‖                               ‖           ||           ‖                               ‖           ||           ‖                               ‖    
        2                               2           ||           2                               2           ||           2                               2           ||           2                               2    

## Excel

The Excel file contains all of the annotated logic gates and circuits. Although each logic function has a separate sheet within the Excel file, there are often multiple versions of each gate (and sometimes versions which use unpaired inputs/outputs). The flow of data through a circuit is mapped using colored lines, and there are comments where necessary to help explain how things work.

## Python

This python script exits to quickly solve Hashiwokakero puzzles, and demonstrate how the Hashiwokakero logic gates I've made work. To use, run the script and then execute commands in the console. The main console commands worth knowing are:
    
gates[channel][name].solve(inputs:str, show:bool, prints:bool) -> str: 

This command solves a specific gate for specific inputs. The inputs and outputs are binary strings (ex. '0110'). When the 'show' control is True, the gate is displayed after it is solved ('show' defaults to False). When the 'prints' control is True, the gate prints whether the puzzle was solved or not, if there was an error, and what the maximum recursion depth of the solver was ('prints' defaults to False). Finally, 'channel' is either 1 or 2, and 'name' is thename of the gate you want to use (valid options are presented below). For example:

    In [x]: gates[2]['or'].solve('01', True)

        2                         2    
        ‖                         ‖    
    2 = 8 ======================= 8 = 2
        ‖                         ‖    
    3 = 7 ------------- 3         ‖    
    |   ‖               ‖         ‖    
    |   6 === 4 = 2     ‖         ‖    
    |   ‖               ‖         ‖    
    3 = 7 - 2 ----- 1   ‖         ‖    
        ‖               2         ‖    
    2 = 8 === 3           2 = 5 = 8 = 2
        ‖     |   2         2 |   ‖    
    3 - 7 = 2 | 1 ‖ 2 = 4   ‖ 3 = 7 - 3
    ‖   ‖     | | ‖     ‖   ‖     ‖   ‖
    ‖   6 === 3 | 4 === 8 = 7 === 6   ‖
    ‖   ‖       |       ‖   |     ‖   ‖
    3 - 7 = 3 - 2       2   3 === 7 - 3
        ‖                         ‖    
    2 = 8 ======================= 8 = 2
        ‖                         ‖    
        2                         2    

    Out[x]: '1'

truthTable(gate, channel:int) -> None:

This command prints a truth table for a given gate. 'gate' can be the name of a gate as a string (ex: 'half adder') in which case you might need to specify the number of channels ('channel' defaults to 2). Alternatively, 'gate' can be a Hashiwokakero puzzle object (usually from the 'gates' dictionary: ex. gates[1]['conditional swap']). For example:

    In [x]: truthTable('buffered nand')

    A B | C
    ----|--
    0 0 | 1
    0 1 | 1
    1 0 | 1
    1 1 | 0


timeGate(gate, N:int, channel:int, Print:bool, Return:bool) :

This command solves the puzzle N times for each unique set of inputs. It then either prints a table of the results, or it returns the raw data (or both). 'gate' and 'channel' serve to fetch gates the same way they do for truthTable(). For example:

    In [x]: timeGate('swap')

    Breakdown by Input:
    00|00:  1.4 ± 0.2 ms
    01|10:  1.2 ± 0.2 ms
    10|01:  0.8 ± 0.1 ms
    11|11:  0.8 ± 0.1 ms

    Average Solve Time: 1.0 ± 0.2 ms

Notably, uncertainties are uncertainties on the mean, not the standard deviation you expect for any given measurement. To compute the measurement-wise standard deviation multiple the uncertainty on the mean by sqrt(N).


analyze(N:int, Return:bool) :

This command times every gate in the 'gates' dictionary (including older versions of gates), and then plots the average solve times by the number of vertices, edges, and bridges in each puzzle. N is fed to timeGate(), and the data is returned when 'Return' is True. For example:

    In [x]: analyze(100)

    [a matplotlib plot opens in a new window]

    Out[x]: np.vstack((vertices, 
                       edges, 
                       bridges, 
                       uncertainty on the number of bridges,
                       mean solve time,
                       unvertainty on the mean solve times))

1-CHANNEL LOGIC GATES:

    not                         : (A) -> (NOT A)
    split                       : (A) -> (A, A, A)
    swap                        : (A, B) -> (B, A)
    and or                      : (A, B) -> (A AND B, A OR B)
    nand nor                    : (A, B) -> (A NAND B, A NOR B)
    or                          : (A, B) -> (A OR B, A OR B, B, A)
    nor                         : (A, B) -> (A NOR B, A NOR B, B, A)
    and                         : (A, B) -> (A AND B, A AND B, B, A)
    nand                        : (A, B) -> (A NAND B, A NAND B, B, A)
    xor v0                      : (A, B) -> (A XOR B) *requires non-local logic
    xor                         : (A, B) -> (A XOR B)
    xnor                        : (A, B) -> (A XNOR B)
    half adder                  : (A, B) -> (A AND B, A AND B, A XOR B)
    conditional swap            : (A, B, C) -> (A, C, B) if A else (A, B, C)

2-CHANNEL LOGIC GATES:

    split                       : (A) -> (A, A)
    double split                : (A) -> (A, A, A)
    not                         : (A) -> (NOT A)
    buffer (_, 2)**             : (A) -> (A)
    swap                        : (A, B) -> (B, A)
    or (_, 2, 3)*               : (A, B) -> (A OR B)
    buffered or**               : (A, B) -> (A OR B)
    nor (_, 2)*                 : (A, B) -> (A NOR B)
    buffered nor**              : (A, B) -> (A NOR B)
    and (_, 2)*                 : (A, B) -> (A AND B)
    buffered and**              : (A, B) -> (A AND B)
    nand (_, 2)*                : (A, B) -> (A NAND B)
    buffered nand**             : (A, B) -> (A NAND B)
    xor                         : (A, B) -> (A XOR B)
    buffered xor**              : (A, B) -> (A XOR B)
    xnor                        : (A, B) -> (A XOR B)
    buffered xnor**             : (A, B) -> (A XNOR B)
    conditional source          : (A, B, C) -> (C) if A else (B)
    conditional source 2        : (B, A, C) -> (C) if A else (B)
    conditional swap            : (A, B, C) -> (C, B) if A else (B, C)
    half adder                  : (X, Y) -> (X AND Y, X XOR Y)
    buffered half adder**       : (X, Y) -> (X AND Y, X XOR Y)
    full adder - bottom***      : (X, Y) -> (Cout, Sum)
    full adder - middle***      : (X, Y, Cin) -> (Cout, Sum)
    full adder - top 1***       : (X, Y, Cin) -> (Cout, Sum)
    full adder - top 2***       : (X, Y, Cin) -> (Sum)
    2b adder****                : (X1, X0, Y1, Y0) -> (S2, S1, S0)
    4b adder****                : (XXXX, YYYY) -> (SSSSS)
    8b adder****                : (XXXXXXXX, YYYYYYYY) -> (SSSSSSSSS)
    compliment - bottom*****    : (X) -> (Cout, X')
    compliment - middle*****    : (X, Cin) -> (Cout, X')
    compliment - top*****       : (X, Cin) -> (X')
    4b compliment*****          : (XXXX) -> (XXXX')
    8b compliment*****          : (XXXXXXXX) -> (XXXXXXXX')
    
    *       Multiple versions of these gates exist. For example, 'or', 'or 2', 
            and 'or 3' are all gates which compute an OR function using slightly 
            different designs.
    **      Buffered gates include a buffer to protect the inputs/outputs from
            any potential setups that could lead to undefined logic. I'm not 
            sure if this is a thing that can actually happen, but the buffered 
            gates are definitely safe. That being said, adding a buffer also 
            adds additional vertices/edges/bridges, which means that they 
            require a little more time to solve.
    ***     There are several versions of full adders, which are each designed 
            to be different parts of a stack of full adders. The bottom full 
            adder doesn't have a carry in, so is actually just a half adder 
            where the carry out lines up with the carry in of the full adder 
            above it. The middle version has aligned carry in's/out's, so you 
            can stack as many of them on top of each other as you need. Finally, 
            there are 2 versions of full adders to go on top of the stack, one 
            which includes a carry out (in a different position than the other 
            full adders), and one which doesn't have a carry out (and is 
            essentially just 2 XOR gates). As these gates are stackable, it is 
            easy to use them (and swap gates) to create n-bit addition circuits 
            like the 4-bit and 8-bit adders (the 2-bit adder is constructed 
            directly from half adders).
    ****    'nb adder's are n-bit adder which computes the sum of two n-bit 
            numbers, which it outputs as an (n + 1)-bit number. I Currently 
            have 2-bit, 4-bit, and 8-bit adders. All of which have been 
            implemented in the function add(x, y, bits) (located at the bottom 
            of the script). This function takes decimal inputs, computes the 
            sum of the inputs using the puzzle specified by the number of bits, 
            and then outputs the number in decimal.
    *****   The 'compliment - ...' gates are pieces which can be stacked together 
            to make an n-bit circuit which computes the 2's compliment of the 
            input number. They work by notting all of the input bits, and then 
            adding 1 to the result (with half adders). The '4b compliment' and 
            '8b compliment' circuits are examples of how the pieces fit together 
            to create a circuit that can compute the compliment of an n-bit 
            number. Also, similar to the add(x, y) function the function 
            subtract(x, y) uses the 8-bit compliment circuit and the 8-bit 
            adder circuits to compute (x - y). Although it it's a single puzzle, 
            it is fun to be able to both add and subtract numbers!

@author: Violet Saathoff
