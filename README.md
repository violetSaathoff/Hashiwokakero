# Hashiwokakero

"""
Created on Sun Nov 28 11:08:53 2021
This script exits to quickly solve Hashiwokakero puzzles, and demonstrate 
how the Hashiwokakero logic gates I've made work. To use, run the script 
and then execute commands in the console. The main console commands worth 
knowing are:
    
    gates[channel][name].solve(inputs:str, show:bool, prints:bool) -> str: 
        This command solves a specific gate for specific inputs. The inputs 
        and outputs are binary strings (ex. '0110'). When the 'show' control 
        is True, the gate is displayed after it is solved ('show' defaults 
        to False). When the 'prints' control is True, the gate prints 
        whether the puzzle was solved or not, if there was an error, and 
        what the maximum recursion depth of the solver was ('prints' defaults 
        to False). Finally, 'channel' is either 1 or 2, and 'name' is the 
        name of the gate you want to use (valid options are presented below).
        For example:
            
            In [x]: gates[2]['or'].solve('01', True)
        
                2                         2    
                ‖                         ‖    
            2 = 8 ======================= 8 = 2
                ‖                         ‖    
            3 = 7 ------------- 3         ‖    
            |   ‖               ‖         ‖    
            F = 8 === 4 = 2     ‖         ‖    
            |   ‖               ‖         ‖    
            3 = 7 - 2 ----- 1   ‖         ‖    
                ‖               2         ‖    
            2 = 8 === 3           2 = 5 = 8 = 2
                ‖     |   2         2 |   ‖    
            3 - 7 = 2 | 1 ‖ 2 = 4   ‖ 3 = 7 - 3
            ‖   ‖     | | ‖     ‖   ‖     ‖   ‖
            T = 8 === 3 | 4 === 8 = 7 === 8 = T
            ‖   ‖       |       ‖   |     ‖   ‖
            3 - 7 = 3 - 2       2   3 === 7 - 3
                ‖                         ‖    
            2 = 8 ======================= 8 = 2
                ‖                         ‖    
                2                         2    
                
            Out[x]: '1'
        
        For proper documentation, run >>> help(Hashiwokakero.solve)
    truthTable(gate, channel:int) -> None:
        This command prints a truth table for a given gate. 'gate' can be 
        the name of a gate as a string (ex: 'half adder') in which case you 
        might need to specify the number of channels ('channel' defaults to 2).
        Alternatively, 'gate' can be a Hashiwokakero puzzle object (usually 
        from the 'gates' dictionary: ex. gates[1]['conditional swap']).
        For example:
            
            In [x]: truthTable('buffered nand')
            
            A B | C
            ----|--
            0 0 | 1
            0 1 | 1
            1 0 | 1
            1 1 | 0
        
        For proper documentation, run >>> help(truthTable)
    
    
    timeGate(gate, N:int, channel:int, Print:bool, Return:bool) :
        This command solves the puzzle N times for each unique set of inputs. 
        It then either prints a table of the results, or it returns the raw 
        data (or both). 'gate' and 'channel' serve to fetch gates the same 
        way they do for truthTable(). For example:
            
            In [x]: timeGate('swap')
            
            Breakdown by Input:
            00|00:  1.4 ± 0.2 ms
            01|10:  1.2 ± 0.2 ms
            10|01:  0.8 ± 0.1 ms
            11|11:  0.8 ± 0.1 ms
            
            Average Solve Time: 1.0 ± 0.2 ms
        
        Notably, uncertainties are uncertainties on the mean, not the standard 
        deviation you expect for any given measurement. To compute the 
        measurement-wise standard deviation multiple the uncertainty on the 
        mean by sqrt(N).
    
    For proper documentation, run >>> help(timeGate)
    
    
    analyze(N:int, Return:bool) :
        This command times every gate in the 'gates' dictionary (including 
        older versions of gates), and then plots the average solve times 
        by the number of vertices, edges, and bridges in each puzzle. N 
        is fed to timeGate(), and the data is returned when 'Return' is True.
        For example:
            
            In [x]: analyze(100)
            
            [a matplotlib plot opens in a new window]
            
            Out[x]: np.vstack((vertices, 
                               edges, 
                               bridges, 
                               uncertainty on the number of bridges,
                               mean solve time,
                               unvertainty on the mean solve times))
        
        For proper documentation, run >>> help(analyze)
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
"""
