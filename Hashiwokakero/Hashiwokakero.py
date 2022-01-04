# -*- coding: utf-8 -*-
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

#Import Libraries
import os
import numpy as np
import time
from collections import OrderedDict
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#Round a Number to a Specific Number of Sig Figs
def sfRound(x, sf):
    """Round a Number to a Certain Number of Significant Figures
    
    Parameters
    ----------
    x : float
        The value to be rounded.
    sf : int
        The number of significant figures to round to.

    Returns
    -------
    x : float
        The rounded value.

    """
    
    #Round x
    x = round(x, int(max(sf - 1, 0) - np.floor(np.log10(abs(x)))))
    
    #Convert x to an Integer (If Applicable)
    if int(x) == x:
        x = int(x)
    
    #Return x
    return x

#Round a Number and its Uncertainty to the Appropriate Precision
def uRound(x, ux):
    """Round a Number and its Uncertainty
    
    For proper reporting, uncertainty should typically be reported to 1 
    significant figure, and the value should be rounded to the same 
    precision as the uncertainty.
    
    Parameters
    ----------
    x : float
        The principle value.
    ux : float
        The uncertainty in x.

    Returns
    -------
    x : float
        The rounded value.
    ux : float
        The rounded uncertainty.

    """
    
    #Get the Number of Decimals to Round to
    n = -int(np.floor(np.log10(abs(ux))))
    
    #Round x and ux
    x, ux = round(x, n), round(abs(ux), n)
    
    #Check if x and ux are Integers
    if ux >= 1:
        x, ux = int(x), int(ux)
    
    #Return x and ux
    return x, ux

#Reduced Chi-Squared Calculator
def chisq(residuals, sigma, ddof = 2):
    """Compute the Reduced Chi-Squared of a Fit

    Parameters
    ----------
    residuals : np.array
        The residuals in y (data - fit line).
    sigma : np.array/float
        The uncertainty in y.
    ddof : int
        The number of delta degrees of freedom (typically equal to the
        number of fit parameters). 
        The default is 2.

    Returns
    -------
    tuple of non-negative floats
        (chi-squared, uncertainty on chi-squared)
    
    Notes
    -----
    Both of the outputs of this function are only valid when there are 
    enough data points that the chi-squared distribution can be well 
    approximated by a guassian distribution (i.e. the more data the 
    better). While this is somewhat reflected in the uncertainty, the 
    actual uncertainty on chi-squared is asymmetric when the number of
    data points is small since the chi-squared distribution itself is
    asymmetric (this is largely because chi-squared >= 0). Since reduced
    chi-squared is always expected to be roughly 1, you can estimate the 
    number of data points you'll need to have n sigma confidence in
    chi-squared using the formula: N = 2*n**2 + ddof where n is how 
    confident you would like to be (in standard deviations), and N is the 
    minimum number of relavent data points you'll need to justify your
    results.
    """
    
    #Compute and Return Chi-Squared
    k = len(residuals[ddof:])
    chi = np.sum((residuals/sigma)**2)/k
    return chi, np.sqrt(2/k)

#A Class for the Verticies
class Vertex:
    """
    A Sub-Class for the Vertices in Hashiwokakero Puzzles
    
    You should not manipulate the Vertices directly.
    
    Attributes
    ----------
    face : int
        The Face Value of the vertex (inputs are given a face 
        value of 9, and outputs are given a face value of 10).
    value : int
        The number of bridges the vertex needs to be connected to 
        for the puzzle to be solved (this number changes as the 
        puzzle is solved)
    index : int
        The index of the vertex withing the Hashiwokakero puzzle's 
        list of vertices.
    position : (int, int)
        The 0-indexed (row, column) position of the vertex within 
        the grid of the Hashiwokakero puzzle.
    neighbors : list(Vertex)
        The other vertices the vertex can build bridges to. This list 
        is paired with the edge list.
    edges : list(Edge)
        The edges which connect to the neighboring vertices. This list 
        is paired with the vertex list.
    
    Methods
    -------
    __init__(self, value:int, index:int, position:(int, int)) -> None : 
        Create a new Vertex Object
    reset(self) -> None : 
        Reset the Vertex to its Initial State
    copy(self, reset:bool = False) -> Vertex : 
        Create a Shallow Copy of the Vertex
    addForced(self, puzzle:Hashiwokakero) -> int : 
        Add bridges which are forced by logical deductions.
    
    """
    
    #Create a New Vertex Object
    def __init__(self, value:int, index:int, position:tuple):
        """Create a New Vertex Object"""
        self.face = value #The Number on the Vertex
        self.value = value #The Number of Bridges Needed to Satisfy the Vertex
        self.index = index #The Index of the Vertex in the Vertices List in the Hashiwokakero Puzzle Object
        self.position = position #The Position of the Vertex Within the Grid
        self.neighbors = [] #The Vertices this Node Connects to (Paired with the Edges)
        self.edges = [] #The Edges which Connect to this Vertex
    
    #Reset the Vertex
    def reset(self):
        """Reset a Vertex to its Initial State"""
        self.value = self.face
    
    #Shallow Copy the Vertix
    def copy(self, reset:bool = False):
        """
        Copy a Vertex
        
        Parameters
        ----------
        reset : bool, optional
            When True, the copy will be reset. When False, the 
            copy will be in the same state as the original vertex. 
            The default is False.

        Returns
        -------
        new : Vertex
            A Copy of the Original Vertex.

        """
        #Make a New Node
        new = Vertex(self.face, self.index, self.position)
        
        #Shallow Copy the Edges
        new.edges = self.edges.copy()
        
        #Shallow Copy the Neighbors
        new.neighbors = self.neighbors.copy()
        
        #Update the Value
        if not reset:
            new.value = self.value
        
        #Return the New Vertex
        return new
    
    #Add Forced Bridges
    def addForced(self, puzzle) -> int:
        """
        Add Forced Bridges Around the Given Vertex
        
        (This is only to be used by the Hashiwokakero solver)
        
        
        Logic Used
        ----------
        
        Saturation : 
            If the number of bridges the vertex needs to be 
            complete is equal to the number of bridges it 
            can build in all directions combined, the every 
            possible bridge must be built.  
            
            For example, a 4 with no bridges and 2 other vertices 
            it can connect to must build double bridges to each of 
            those vertices.
        
        Semi-Saturation : 
            If the number of bridges the vertex needs to be 
            complete is one less than the total number of 
            bridges it can build in all directions combined, 
            then every edge where a double bridge could be built 
            must at least have a single bridge.
            
            For example, a 5 with not bridges and 3 other vertices 
            it can connect to build build a single bridge in each 
            direction (2 of those edged will later become double 
            bridges, but we don't know which 2 at this point).
        
        Over-Saturation : 
            If the number of bridges the vertex needs to be 
            complete exceedes the total number of bridges it 
            can build in all directions combined, then the 
            current vertex is not solvable (implying that 
            either the puzzle is not solvable, or the solver 
            made an incorrect choice before it recursed).
        
        Puzzle Connectivity : 
            Unless there are only 2 vertices in the puzzle, a 
            1 cannot build a bridge to another 1, and a 2 cannout 
            build a double bridge to another 2, since otherwise 
            the finaly puzzle would not be fully connected (this 
            is actually implemented before the other logical tools 
            detailed since it makes it more likely that those 
            tools will be relavent).
        
        """
        
        #Get How Many Bridges Can be Added in Each Direction
        possible = [edge.canBuild() for edge in self.edges]
        
        #Count the Number of Bridges that Could be Build
        N = sum(possible)
        M = sum(x > 0 for x in possible)
        
        #Remove Forbidden Bridges
        if self.face == 1:
            #Look for Other 1's
            for i, vertex in enumerate(self.neighbors):
                if vertex.face == 1 == possible[i]:
                    #A Bridge Between 1's is Forbidden
                    possible[i] = 0
                    N -= 1
                    M -= 1
        elif self.face == 2:
            #Look for Other 2's
            for i, vertex in enumerate(self.neighbors):
                if vertex.face == 2 == possible[i]:
                    #A Double Bridge Between 2's is Forbidden
                    possible[i] = 1
                    N -= 1
        
        #Check the Case
        if self.face == 10:
            #OUTPUT: Write the Output if Possible
            if N == 0:
                #Count the Number of Bridges Connected to the Output
                count = 10 - self.value
                
                #Update the Puzzle's Sum
                puzzle.sum += count
                
                #Record the Truth Value of the Output
                self.face = [0,'F','T',3,'F',5,'T',7,8][count]
                
                #Mark that it Requires No More Bridges
                self.value = 0
                
                #Record the Move
                puzzle.moves.append((1, self, count))
            
            #Return that an Output was Marked
            return 2
        elif self.value > N:
            #Check for an Edge Case
            if len(puzzle.vertices) == 2 and self.edges and 3 > self.face == self.neighbors[0].face:
                """
                The edge case in question is when the puzzle consists of 
                2 verticies which are both either 1's or 2's. In those 
                cases my solver will initially assume that it is forbidden 
                to add the bridges required to solve the puzzle, and will 
                then see that it can't add enough bridges to solve the 
                puzzle. To get around this I wait until there is an 
                error (a rare enough occurence that this rarely affects 
                the run time at all), then I check to see if the puzzle 
                is in the edge case. And if it is I add the bridges required 
                to solve the puzzle.
                """
                #Get the Number of Bridges
                bridges = self.edges[0].canBuild()
                
                #Add Bridges
                self.edges[0].build(bridges)
                
                #Record the Move
                puzzle.moves.append((0, self.edges[0], bridges))
                
                #Return that Bridges Were Made
                return 1
            
            #OVER SATURATION: No Solution is Possible
            puzzle.moves.append((-1, self))
            
            #Return an Error Code
            return -1
        elif N == self.value:
            #SATURATION: Add All Possible Bridges
            for edge, bridges in zip(self.edges, possible):
                #Only Add Bridges Which are Possible to Add
                if bridges:
                    #Build the Bridges
                    edge.build(bridges)
                    
                    #Record the Move
                    puzzle.moves.append((0, edge, bridges))
            
            #Return the Bridges Were Added
            return 1
        elif self.value == N - 1 and N > M:
            #SEMI-SATURATION: Add 1 Bridge in Each Direction With Degeneracy
            for edge, bridges in zip(self.edges, possible):
                if bridges == 2:
                    #Build the Bridges
                    edge.build(1)
                    
                    #Record the Move
                    puzzle.moves.append((0, edge, 1))
            
            #Return the Bridges Were Added
            return 1
        
        #Return that No Bridges Were Made
        return 0

#A Class for the Edges
class Edge:
    """
    A Sub-Class for the Edges in a Hashiwokakero Puzzle
    
    You should not manipulate the Edges directly.
    
    Attributes
    ----------
    A : Vertex
        The first vertex (in reading order) which the edge connects to.
    B : Vertex
        The second vertex (in reading order) which the edge connects to.
    index : int
        The index of the edge in the puzzle's list of edges.
    orientation : bool
        Whether or not the edge is horizontal.
    bridges : int
        When non-negative, 'bridges' is the number of bridges currently 
        built along the edge. When negative, abs(self.bridges) is the 
        number of other edges which are currently blocking the edge.
    intersections : list(Edge) 
        A list of Edges which intersect with this edge. When an edge 
        builds its first bridge it updates the bridge count of all 
        edges which intersect with it to reflect the fact that those 
        edges are now blocked and cannot build bridges. When an edge 
        destroys its last bridge it updates the bridge coutn of all 
        edges which intersect with it to reflect the fact that those 
        edges are no longer blocked, and can build bridges again 
        (unless they're blocked by a different edge). This system 
        means that determining if an edge is blocked is an O(1) 
        operation, but that building/destroying a bridge might run 
        in O(intersections) time, or it might run in O(1) time. This 
        is nice though, because checking if a connection is open is 
        the far more common method to use (every time the puzzle 
        investigates a vertex, it runs edge.canBuild() for each edge 
        connected to that vertex, and then maybe it builds bridges 
        along some of the edges).
    fixed : bool
        When True, no more bridges are allowed to be built along an 
        edge even if they could be. This allows inputs to be set 
        in 2-channel puzzles that don't use blank nodes.
    
    Methods
    -------
    __init__(self, A:Vertex, B:Vertex, index:int) -> None : 
        Create a New Edge Object
    reset(self) -> None : 
        Reset the Edge to its Initial State
    copy(self, reset:bool = False) -> Edge : 
        Create and Return a Shallow Copy of the Edge
    canBuild(self) -> int : 
        Return the Number of Bridges Which can be Built Along the Edge
    build(self, bridges:int) -> None : 
        Build the Specified Number of Bridges Along the Edge
    destroy(self, bridges:int) -> None : 
        Destroy the Specified Number of Bridges Along the Edge
    
    """
    
    #Create a New Edge Object
    def __init__(self, A:Vertex, B:Vertex, index:int):
        """Create a New Edge Object"""
        self.A = A #The Vertex which is more Up/Left
        self.B = B #The Vertex which is more Down/Right
        self.index = index #The Index of the Edge in the Edges List in the Hashiwokakero Puzzle Object
        self.orientation = A.position[0] == B.position[0] #Whether the Edge is Horizontal or Not (Vertical)
        self.bridges = 0 #The Number of Bridges Currently Along the Edge (Or when Negative, the Number of Edges Blocking this Edge)
        self.intersections = [] #Edges which Intersect with this Edge
        self.fixed = False #When True, Additional Bridges are Not Allowed to be Built
    
    #Reset the Edge
    def reset(self) -> None:
        """Reset the Edge to its Initial State"""
        self.bridges = 0
        self.blocks = 0
        self.fixed = False
    
    #Shallow Copy the Edge
    def copy(self, reset:bool = False):
        """
        Copy an Edge
        
        Parameters
        ----------
        reset : bool, optional
            When True, the copy will be reset. When False, the 
            copy will be in the same state as the original edge. 
            The default is False.
        
        Returns
        -------
        new : Edge
            A Copy of the Original Edge.

        """
        
        #Make a New Node
        new = Edge(self.A, self.B, self.index)
        
        #Copy the Intersection Indicies (This can only be updated after all edges are copied)
        new.intersections = self.intersections.copy()
                
        #Update the State of the New Vertex
        if not reset:
            new.value = self.bridges
        
        #Return the New Edge
        return new
    
    #Determine How Many Bridges Can be Built - O(1)
    def canBuild(self) -> int:
        """Compute How Many Bridges Can Currently be Built Along this Edge"""
        return 0 if self.bridges < 0 or self.fixed else min((self.A.value, self.B.value, 2 - self.bridges))
    
    #Build Bridges - O(intersections)/O(1)
    def build(self, bridges:int) -> None:
        """
        Build the Specified Number of Bridges Along this Edge
        
        Parameters
        ----------
        bridges : int
            The number of bridges to be built.
        
        Notes
        -----
        There are no safety checks to make sure the bridges can 
        be built, so this is really only for the solver to use.
        """
        
        #Mark Intersecting Edges as Blocked
        if not self.bridges:
            for edge in self.intersections:
                edge.bridges -= 1
        
        #Update the Number of Bridges Along this Edge
        self.bridges += bridges
        
        #Update the Values of the Vertices at Each End
        self.A.value -= bridges
        self.B.value -= bridges
    
    #Destroy Bridges
    def destroy(self, bridges:int) -> None:
        """
        Destroy the Specified Number of Bridges Along this Edge
        
        Parameters
        ----------
        bridges : int
            The number of bridges to be destroyed.
        
        Notes
        -----
        There are no safety checks to make sure the bridges can 
        be destroyed, so this is really only for the solver to use.
        """
        
        #Update the Number of Bridges Along this Edge
        self.bridges -= bridges
        
        #Update the Values of the Vertices at Each End
        self.A.value += bridges
        self.B.value += bridges
        
        #Mark Intersecting Edges as Unblocked
        if not self.bridges:
            for edge in self.intersections:
                edge.bridges += 1

#A Class for the Puzzle
class Hashiwokakero:
    """
    A Hashiwokakero Puzzle Class
    
    Attributes
    ----------
    grid : np.ndarray
        The input grid which was used to initialize the puzzle. This grid 
        is never referenced after initialization, but is stored anyway.
    shape : (int, int)
        The shape of the puzzle (rows, columns).
    size : (int, int, int)
        Measures of the size of the puzzle besides the physical size of 
        the grid: (#vertices, #edges, #bridges to solve)
    type : int
        What kind of logic the puzzle uses:
            0 : None   - the puzzle isn't a circuit)
            1 : (1, 2) - the puzzle is a 1-channel circuit
            2 : (4, 6) - the puzzle is a 2-channel circuit
    logic : None or (int, int)
        The logic levels associated with the puzzle's type (low, high).
    vertices : list(Vertex)
        A list of the Vertices in the puzzle.
    edges : list(Edge)
        A list of the Edges in the puzzle.
    queue : OrderedDict(int)
        A priority queue of active Vertex indices (an active vertex is 
        one which isn't full).
    backup : list
        A record of how indices were removed from the queue which 
        facilitates backtracking
    moves : list
        A record of how bridges were added to the puzzle which 
        facilitates backtracking.
    inputs : list(Vertex)
        A list of input Vertices
    outputs : list(Vertex)
        A list of output Vertices
    sum : int
        Twice the number of bridges needed to solve the puzzle.
    SUM : int
        The initial value of the sum.
    solved : bool
        A flag which markes if the puzzle was sucessfully solved.
    error : bool
        A flag which marks if the puzzle is in an unsolvable state.
    
    
    Methods
    -------
    __init__(self, grid:np.ndarray) -> None:
        Initialize a new Hashiwokakero puzzle object from a 
        2D grid representation of the puzzle.
    __repr__(self) -> None : 
        Creates an ASCII representation of the puzzle in its current state.
    reset(self) -> None : 
        Resets the puzzle to its initial state.
    copy(self) -> Hashiwokakero :
        Deep copies the puzzle.
    writeInputs(self, inputs:str) -> None : 
        Writes the given inputs to the inputs of the puzzle.
    readOutputs(self) -> str : 
        Reads the outputs of the puzzle.
    undo(self) -> None : 
        Removes the last bridge which was created.
    connected(self) -> bool : 
        Checks if the puzzle is fully connected.
    solveForced(self): -> bool : 
        Solves all logically-forced bridges from the current state.
    key(self, index:int) -> int : 
        A value function for the recursive solver.
    solver(self, depth:int = 0, maxDepth:int = 0) -> int : 
        A recursive solver which manages how the puzzle is solved.
    solve(self, inputs:str = '', show:bool = False, prints:bool = False) -> str : 
        A method which writes the puzzle inputs (using self.writeInputs()), 
        solves the puzzle (using self.solver()), reads the outputs (using 
        self.readOutputs()), and then determines what information to print to 
        the console based on the 'show'/'prints' controls.
    save(self, filename:str) -> None : 
        Saves the puzzle to the given filename/filepath as a .csv file.
    render(self, reset:bool = False) -> None : 
        Creates a JPEG rendering of the puzzle for professional display.
    
    """
    
    #Integer Type Pool
    integers = {int, np.int8, np.int16, np.int32, np.int64}
    strings = {str, np.str_}
    
    #Create a New Puzzle Object
    def __init__(self, grid:np.ndarray) -> None:
        """Make a New Hashiwokakero Object"""
        #Initialize the Attributes
        self.grid = grid #The Puzzle Grid (saved for convenience, this is never actually used)
        self.shape = self.grid.shape #The Shape of the Puzzle
        self.size = None
        self.type = 0 #What Sorts of Inputs/Outputs the Puzzle Uses
        self.logic = None #The Logic Levels Associated with the Type of Inputs/Outputs
        self.vertices = [] #A List of the Vertices in the Puzzle
        self.edges = [] #A List of the Edges which Connect the Vertices
        self.queue = OrderedDict() #A Queue of Vertex Indices for Solving the Puzzle
        self.backup = [] #A Backup Record for the Queue (Used when Backtracking)
        self.moves = [] #A List of Which Bridges were Added in Which Order
        self.inputs = [] #A List of Input Vertices
        self.outputs = [] #A List of Output Vertices
        self.sum = 0 #The Number of Faces of the Vertices (Twice the Number of Bridges Required to Solve the Puzzle)
        self.SUM = 0 #A Backup for self.sum (used when resetting the puzzle)
        self.solved = False #A Flag for if the Puzzle is Solved
        self.error = False #An Error Flag
        
        #Find the Vertices While Adding Edges
        columns = [None]*self.shape[1]
        edgemap = -np.ones(self.shape, int)
        vertical = -np.ones(self.shape, int)
        for i, row in enumerate(grid):
            last = None
            for j, square in enumerate(row):
                #Check for a Possible New Vertex
                if type(square) in self.integers or (type(square) in self.strings and square.isdigit()):
                    #Cast the Square to an Integer
                    square = int(square)
                    
                    #Check if it's Actually a Vertex
                    if 0 < square < 11:
                        #Make the New Vertex
                        vertex = Vertex(square, len(self.vertices), (i, j))
                        
                        #Record the Vertex
                        self.vertices.append(vertex)
                        
                        #Add the Vertex to the Queue
                        self.queue[vertex.index] = -1
                        self.backup.append(-1)
                        if 7 <= vertex.value <= 8:
                            #Move it to the Front
                            self.queue.move_to_end(vertex.index, False)
                        
                        #Record if it's an Input or an Output
                        if vertex.value == 9:
                            #Record the Input
                            self.inputs.append(vertex)
                        elif vertex.value == 10:
                            #Record the Output
                            self.outputs.append(vertex)
                        else:
                            #Update the Sum
                            self.sum += vertex.value
                        
                        #Add a Horizontal Edge to the Previous Vertex
                        if last:
                            #Make the Edge
                            edge = Edge(last, vertex, len(self.edges))
                            
                            #Record the Edge
                            self.edges.append(edge)
                            
                            #Add the Edge to the Vertices
                            last.edges.append(edge)
                            last.neighbors.append(vertex)
                            vertex.edges.append(edge)
                            vertex.neighbors.append(last)
                            
                            #Mark the Edge in the Edge Map
                            for k in range(last.position[1] + 1, j):
                                edgemap[i][k] = edge.index
                        
                        #Add a Vertical Edge to the Previous Vertex
                        if columns[j]:
                            #Make the Edge
                            edge = Edge(columns[j], vertex, len(self.edges))
                            
                            #Record the Edge
                            self.edges.append(edge)
                            
                            #Add the Edge to the Vertices
                            columns[j].edges.append(edge)
                            columns[j].neighbors.append(vertex)
                            vertex.edges.append(edge)
                            vertex.neighbors.append(columns[j])
                            
                            #Find Intersections
                            for k in range(columns[j].position[0] + 1, i):
                                #Mark the Vertical Edge Map
                                vertical[k][j] = edge.index
                                
                                #Check for an Intersection
                                if edgemap[k][j] > -1:
                                    #Record the Intersection
                                    self.edges[edgemap[k][j]].intersections.append(edge)
                                    edge.intersections.append(self.edges[edgemap[k][j]])
                        
                        #Save the Vertex
                        last = vertex
                        columns[j] = vertex
                    elif square == 11:
                        #Update the Puzzle's Type
                        self.type = 3
                        
                        #Save the Position of the Input
                        self.inputs.append((i, j))
                    elif square == 12:
                        #Update the Puzzle's Type
                        self.type = 3
                        
                        #Save the Position of the Output
                        self.outputs.append((i, j))
        
        #Save the Sum
        self.SUM = self.sum
        
        #Determine the Input/Output Type
        if self.inputs:
            #Check the Type
            if self.type == 3:
                #Replace the Inputs With Their Edges
                for index, (i, j) in enumerate(self.inputs):
                    #Check if it's Vertical or Horizontal
                    if vertical[i, j] > -1:
                        #Check for a Conflict
                        if edgemap[i, j] > -1:
                            #Raise a Value Error
                            raise ValueError('Conflicting input edges at ' + str((i, j)))
                        else:
                            #Update the Input
                            self.inputs[index] = self.edges[vertical[i, j]]
                    elif edgemap[i, j] > -1:
                        #Update the Input
                        self.inputs[index] = self.edges[edgemap[i, j]]
                    else:
                        #Raise a Value Error
                        raise ValueError('Disconnected input at ' + str((i, j)))
                
                #Replace the Outputs With Their Edges
                for index, (i, j) in enumerate(self.outputs):
                    #Check if it's Vertical or Horizontal
                    if vertical[i, j] > -1:
                        #Check for a Conflict
                        if edgemap[i, j] > -1:
                            #Raise a Value Error
                            raise ValueError('Conflicting output edges at ' + str((i, j)))
                        else:
                            #Update the Output
                            self.outputs[index] = self.edges[vertical[i, j]]
                    elif edgemap[i, j] > -1:
                        #Update the Output
                        self.outputs[index] = self.edges[edgemap[i, j]]
                    else:
                        #Raise a Value Error
                        raise ValueError('Disconnected output at ' + str((i, j)))
                
                #Save the Logic Levels
                self.logic = (1, 2)
            else:
                #Look at Each Input/Output to Check for Consistency
                types = []
                for vertex in self.inputs + self.outputs:
                    #Count the Number of Edges Which Lead to 3's
                    count = 0
                    for edge in vertex.edges:
                        count += edge.A.face == 3 or edge.B.face == 3
                    
                    #Check the Count
                    if count == 0:
                        types.append(1)
                    elif count == 2:
                        types.append(2)
                    else:
                        raise ValueError('Unrecognized input style at: %s' % str(vertex.position))
                
                #Check for Consistency
                if all(x == types[0] for x in types):
                    #Set the Type
                    self.type = 1 if types[0] == 1 else 2
                    
                    #Set the Logic Levels
                    self.logic = [(1, 2), (4, 6)][self.type - 1]
                else:
                    #Flad the Inconsistency
                    raise ValueError('Inconsistent input styles detected')
        
        #Save the Size of the Puzzle
        b = int(round((self.sum + (len(self.inputs) + len(self.outputs))*[0, 1.5, 5, 0][self.type])/2))
        self.size = (len(self.vertices), len(self.edges), b)
    
    #Print the Puzzle
    def __repr__(self) -> None:
        """Create an ASCII Representation of the Puzzle"""
        #Make a Grid
        grid = np.full(self.shape, ' ', str)
        
        #Fill the Vertices
        for vertex in self.vertices:
            grid[vertex.position] = ['I','O'][vertex.face - 9] if vertex.face in (9, 10) else str(vertex.face)
        
        #Fill the Edges
        for edge in self.edges:
            #Check if There Are Bridges
            if edge.bridges > 0:
                #Check the Orientation
                if edge.orientation:
                    #Add the Horizontal Bridge
                    i, J = edge.B.position
                    char = ['-','='][edge.bridges - 1]
                    for j in range(edge.A.position[1] + 1, J):
                        grid[i][j] = char
                else:
                    #Add the Vertical Bridge
                    I, j = edge.B.position
                    char = ['|','‖'][edge.bridges - 1]
                    for i in range(edge.A.position[0] + 1, I):
                        grid[i][j] = char
        
        #Horizontal Bridges
        horizontal = {'-','='}
        
        #Get 2 Characters at a Time
        def chars(a:str, row:np.ndarray, j:int):
            #Check the Case
            if j == len(row):
                return a
            elif a in horizontal and a == row[j]:
                return 2*a
            else:
                return a + ' '
        
        #Build the String
        return '\n'.join(
            ''.join(chars(char, row, j) for j, char in enumerate(row, start = 1))
            for row in grid)
    
    #Reset the Puzzle
    def reset(self):
        """Reset the Puzzle to its Original State"""
        #Reset the Vertices
        for vertex in self.vertices:
            vertex.reset()
        
        #Reset the Edges
        for edge in self.edges:
            edge.reset()
        
        #Reset the Inputs
        for vertex in self.inputs:
            vertex.face = 9
            vertex.value = 9
        
        #Reset the Outputs
        for vertex in self.outputs:
            vertex.face = 10
            vertex.value = 10
        
        #Reset the Queue
        self.queue = OrderedDict()
        for i in range(len(self.vertices)):
            #Record the Index
            self.queue[i] = -1
            self.backup[i] = -1
            
            #Check the Value
            if 7 <= vertex.value <= 8:
                #Move it to the Front
                self.queue.move_to_end(vertex.index, False)
        
        #Reset the Sum
        self.sum = self.SUM
        
        #Reset the Moves List
        self.moves.clear()
        
        #Reset the Flags
        self.solved = False
        self.error = False
    
    #Copy the Puzzle
    def copy(self, reset:bool = False):
        """Deep Copy the Puzzle"""
        #Make a New Puzzle
        new = Hashiwokakero(self.grid)
        
        #Copy the State
        if reset:
            #Copy the Vertices
            for vertex in self.vertices:
                new.vertices[vertex.index].face = vertex.face
                new.vertices[vertex.index].value = vertex.value
            
            #Copy the Edges
            for edge in self.edges:
                new.edges[edge.index].bridges = edge.bridges
            
            #Copy Easy Attributes
            new.queue = self.queue.copy()
            new.backup = self.backup.copy()
            new.sum = self.sum
            new.solved = self.solved
            new.error = self.error
            
            #Copy the Moves
            for code, *move in self.moves:
                #Check the Code
                if code == 0:
                    new.moves.append((0, new.edges[move[0].index], move[1]))
                elif code == 1:
                    new.moves.append((1, new.vertices[move[0].index], move[1]))
                elif code == -1:
                    new.moves.append((-1, new.vertices[move[0].index]))
        
        #Return the Copy
        return new
    
    #Write the Inputs of a Puzzle
    def writeInputs(self, inputs:str):
        """
        Write the Puzzle's Inputs
        
        Parameters
        ----------
        inputs : str
            The inputs as a binary string (in reading order).
        
        """
        
        #Check the Number of Inputs
        if len(inputs) != len(self.inputs):
            raise ValueError('this puzzle expected %d inputs. %d inputs were provided.' % (len(self.inputs), len(inputs)))
        
        #Check the Type of the Puzzle
        if 0 < self.type < 3:
            #Write the Inputs
            for vertex, value in zip(self.inputs, map(int, inputs)):
                #Mark the Truth Value
                vertex.face = ['F','T'][value]
                
                #Mark the Number of Bridges Needed
                vertex.value = self.logic[value]
        elif self.type == 3:
            #Write the Inputs
            for edge, value in zip(self.inputs, map(int, inputs)):
                #Build the Appropriate Number of Bridges
                edge.build(self.logic[value])
                
                #Record the Move
                self.moves.append((0, edge, edge.value))
                
                #Mark the Edge as an Active Input
                edge.fixed = True
    
    #Read the Outputs of a Puzzle
    def readOutputs(self) -> str:
        """
        Read the Puzzle's Outputs
        
        Returns
        -------
        str
            The outputs as a binary string (in reading order).
        
        """
        
        #Check the Type of the Puzzle
        if self.type == 0:
            #There Are No Outputs
            return None
        elif self.type < 3:
            #Get/Return the Outputs
            return ''.join('1' if vertex.face == 'T' else ('0' if vertex.face == 'F' else 'X') for vertex in self.outputs)
        elif self.type == 3:
            #Get/Return the Outputs
            return ''.join(str(edge.bridges - 1) for edge in self.outputs)
    
    #Undo the Last Bridge Played
    def undo(self):
        """Destroy the Last Bridge that was Created"""
        #Check the Moves List
        if self.moves:
            #Get the Last Move
            code, *move = self.moves.pop()
            
            #Clear the Flags
            self.solved = False
            self.error = False
            
            #Check the Code
            if code == 0:
                #Split the Input
                edge, bridges = move
                
                #Add Nodes Back to the Queue
                for i in (edge.A.index, edge.B.index):
                    self.queue[i] = self.backup[i]
                
                #Destroy the Bridges
                edge.destroy(bridges)
            elif code == 1:
                #Split the Input
                vertex, count = move
                
                #Reset the Output
                vertex.face = 10
                vertex.value = 10 - count
    
    #Check if the Puzzle is Connected
    def connected(self) -> bool:
        """Check if the Puzzle is Fully Connected"""
        
        #Try to Visit All the Vertices
        queue = {0}
        seen = set()
        while queue:
            #Get the Next Vertex From the Queue
            i = queue.pop()
            
            #Mark it as Seen
            seen.add(i)
            
            #Visit its Connected Neighbors
            for edge, neighbor in zip(self.vertices[i].edges, self.vertices[i].neighbors):
                #Check if the Edge has a Bridge and that the Neighbor is New
                if edge.bridges > 0 and neighbor.index not in seen:
                    #Add the Neighbor to the Queue
                    queue.add(neighbor.index)
        
        #Check if all the Vertices were Visited
        return len(seen) == len(self.vertices)
    
    #Solve Forced Bridges
    def solveForced(self) -> bool:
        """Add All Bridges which are Forced by Logical Deductions"""
        #Continue Until No More Bridges Can be Added
        while self.queue:
            #Get the Next Vertex From the Queue
            i, l = self.queue.popitem(False)
            
            #Check if No Progress Was Made Since the Last Time it was Visited
            if l == len(self.moves):
                #Add the Index Back Onto the Queue
                self.queue[i] = l
                
                #Return that the Puzzle Has Not Been Solved
                return False
            
            #Save the Current Number of Moves Made
            l = len(self.moves)
            
            #Add Forced Bridges Around the Current Vertex
            code = self.vertices[i].addForced(self)
            
            #Check the Result
            if code == -1:
                #Return that the Puzzle Has Not Been Solved
                return False
            elif code == 1:
                #Move Interesting Vertices to the Front of the Queue
                for code, edge, bridges in self.moves[l:]:
                    #Move the Intersecting Edges
                    for intersecting in edge.intersections:
                        for vertex in (intersecting.A, intersecting.B):
                            if vertex in self.queue:
                                self.queue.move_to_end(vertex.index, False)
                    
                    #Move the Nodes Connected to the Current Edge
                    for vertex in (edge.A, edge.B):
                        if vertex in self.queue:
                            #Check if the Vertex Should be Removed
                            if vertex.value:
                                #Move it to the Front
                                self.queue.move_to_end(vertex.index, False)
                            else:
                                #Remove it
                                self.queue.pop(vertex.index)
            
            #Add the Current Vertex to the End of the Queue (if it's still active)
            if self.vertices[i].value:
                #Record the Current Number of Moves Made
                self.queue[i] = len(self.moves)
                self.backup[i] = len(self.moves)
        
        #Check if the Puzzle Was Solved
        self.solved = self.connected()
        
        #Return the Solve State
        return self.solved
    
    #A Value Function for the Recursive Solver
    def key(self, index:int):
        """A Value Function for the Recursive Solver"""
        #Get the Vertex
        vertex = self.vertices[index]
        
        #Compute/Return the Cost Function
        return (vertex.value, sum(bool(edge.canBuild()) for edge in vertex.edges))
    
    #A Recursive Puzzle Solver
    def solver(self, depth:int = 0, maxDepth:int = 0) -> int:
        """A Recursive Solver used by the .solve() Method"""
        #Update the Max Depth
        maxDepth = max(depth, maxDepth)
        
        #See if Solving Forced Bridges Solves the Puzzle
        if self.solveForced() or self.error or not self.queue:
            #Return the Maximum Recursion Depth
            return maxDepth
        else:
            #Find a Good Vertex to Test
            vertex = self.vertices[min(self.queue.keys(), key = self.key)]
            
            #Record the Current Number of Moves
            l = len(self.moves)
            
            #Try Adding a Bridge in Each Direction
            for edge in vertex.edges:
                #Check if the Edge is Open
                if edge.canBuild():
                    #Add a Bridge
                    edge.build(1)
                    
                    #Record the Move
                    self.moves.append((0, edge, 1))
                    
                    #Recurse
                    maxDepth = self.solver(depth + 1, maxDepth)
                    
                    #Check if the Puzzle was Solved
                    if not self.queue and self.connected():
                        #Mark the Puzzle as Solved
                        self.solved = True
                        
                        #Return the Maximum Depth
                        return maxDepth
                    
                    #Backtrack
                    for i in range(l, len(self.moves)):
                        self.undo()
            
            #Return the Maximum Recursion Depth
            return maxDepth
    
    #Solve the Puzzle
    def solve(self, inputs:str = '', show:bool = False, prints:bool = False) -> str:
        """
        Solve the Hashiwokakero Puzzle/Circuit
        
        Parameters
        ----------
        inputs : str, optional
            The inputs you wish to provide the circuit. The number 
            of inputs must match the number of inputs for the 
            circuit in question, and the inputs should be provided 
            as a binary string (in reading order with respect to 
            the physical layout of the puzzle).
            The default is '' (No Inputs).
        show : bool, optional
            When True, an ascii rendering of the solved puzzle 
            will be displayed in the console. 
            The default is False.
        prints : bool, optional
            When True, a status update will be printed after the 
            solve detailing if the puzzle was solved, and what 
            the maximum recursion depth of the solver was. 
            The default is False.
        
        Returns
        -------
        str
            The Puzzle's Outputs (in reading order). Like with 
            the inputs, the outputs are a binary string.
        
        """
        
        #Reset the Puzzle
        self.reset()
        
        #Check for Parity Errors
        if self.type != 1 and self.sum%2:
            return 'parity error'
        
        #Write the Inputs
        self.writeInputs(inputs)
        
        #Use the Recursive Solver
        depth = self.solver()
        
        #Read the Outputs
        outputs = self.readOutputs() if not self.error else 'error'
        
        #Print the Result
        if prints:
            if self.solved:
                print('Puzzle Solved! Maximum Recursion Depth: %d' % depth)
            elif self.error:
                vertex = self.moves[-1][1]
                print('Error: the %s at %s is oversaturated.' % (str(vertex.face, vertex.position)))
            else:
                print('Puzzle not Solved.')
        
        #Display the Puzzle
        if show:
            print('')
            print(self, end = '\n\n')
        
        #Return the Outputs
        return outputs
    
    #Save a Gate as a .csv File
    def save(self, filename:str) -> None:
        '''
        Save the puzzle as a csv file.
        
        Parameters
        ----------
        filename : str
            The name you want to save the file as. You can include a 
            file path (assuming it's a valid destination), but you 
            shouldn't include a suffix (all files are saved as .csv 
            files, so anything after the first period you include 
            in your filename will be ignored).
        
        '''
        
        #Make the Grid
        grid = [['']*self.shape[1] for i in range(self.shape[0])]
        for vertex in self.vertices:
            grid[vertex.position[0]][vertex.position[1]] = str(vertex.face)
        
        #Make Sure the Inputs are Correct
        for vertex in self.inputs:
            grid[vertex.position[0]][vertex.position[1]] = '9'
        
        #Make Sure the Outputs are Correct
        for vertex in self.outputs:
            grid[vertex.position[0]][vertex.position[1]] = '10'
        
        #Open the File
        file = open(os.getcwd() + '\\' + filename.split('.')[0] + '.csv', 'w')
        
        #Write the File
        file.write('\n'.join(','.join(row) for row in grid))
        
        #Close the File
        file.close()
    
    #Render as Puzzle as a JPEG in its Current State
    def render(self, reset:bool = False):
        """
        #Render the Puzzle as a JPEG
        
        Parameters
        ----------
        reset : bool, optional
            When True, the puzzle will be 
            reset before it is rendered. 
            The default is False.
        
        Returns
        -------
        None.

        """
        raise NotImplementedError()

#Load a Gate from a Filename
def load(name):
    """Load a Hashiwokakero Logic Gate by Name"""
    
    #Split the Name
    name = name.split('\\')
    
    #Get the Path/Name
    path = '\\'.join(name[:-1])
    name = '.'.join([name[-1].split('.')[0], 'csv'])
    
    #Check if a Path was Specified
    if len(path) == 0:
        #Try All Folders
        for folder in ['2-Channel','1-Channel']:
            if name in os.listdir(folder):
                #Make/Return the Puzzle
                return Hashiwokakero(np.loadtxt('\\'.join([folder, name]), str, None, ','))
    else:
        #Make Sure the File is Where it Says it is
        if name in os.listdir(path):
            #Make/Return the Puzzle
            return Hashiwokakero(np.loadtxt('\\'.join([path, name]), str, None, ','))
    
    #Raise an Exception
    if len(path) > 0:
        path += '\\'
    raise OSError(''.join([path, name, ' not found.']))

#Recursively Search a Folder Tree
def search(target, case:bool = False, root:str = None, hits:None = None):
    """Search a Folder Tree for a Target
    
    Parameters
    ----------
    target : str
        The item you're searching for.
    case : bool, optional
        Whether or not the search should be case sensitive. 
        The default is False.
    root : str, optional
        The root folder. 
        The default is None (current working directory).
    hits : None, optional
        This is necessary to make the search work, I recommend that you
        don't tamper with it. 
        The default is None.
    
    Returns
    -------
    hits : list
        A list of all the found matches.
    
    Notes
    -----
    Setting 'hits' to anything other than a list will have 'hits' immediately
    replaced by an empty list. This new empty list will then be filled with 
    the matches this function finds. If you instead pass a list for 'hits', 
    the matches this function finds will be appended to the list you passed.
    This typically won't break anything, but you could end up having to dig 
    through a much larger list than would otherwise be necessary.
    
    """
    
    #Inialize the Hits
    if type(hits) != list:
        hits = []
    
    #Get the Directory
    if root is None:
        directory = os.listdir()
        path = ''
    elif type(root) == str:
        directory = os.listdir(root)
        path = root
    
    #Search the Folder
    for file in directory:
        #Make the Filename
        if len(path) > 0:
            filename = '\\'.join([path, file])
        else:
            filename = file
        
        #Check the File
        if (case == False and target.lower() in file.lower()) or target in file:
            hits.append(filename)
        
        #Check if the File is a Folder
        if '.' not in file:
            #Recurse
            hits = search(target, case, filename, hits)
    
    #Return the Hits
    return hits
    
#Load the Gates
def loadGates():
    """Load Hashiwokakero Logic Gates

    Returns
    -------
    gates : dict
        A dictionary of all of the gates in the current working directory.
        The dictionary is organized identically to the folder tree it 
        searched, except that folders names 'n-Channel' where n is an integer
        are indexed using the integer n. For example, gates[1]['xor'] would be
        an xor gate found in the folder '1-Channel', and gates[2]['OLD']['or']
        would be a gate found in the 'OLD' subfolder within the folder 
        '2-Channel'.

    """
    
    #Search the Folder Tree
    hits = search('.csv', False, None, None)
    
    #Initialize the Dictionary
    gates = {}
    
    #Load the Files
    for filename in hits:
        #Split the Filename
        split = filename.split('\\')
        
        #Try Block for Safety
        try:
            #Initialize Sub-Dictionaries
            dictionary = gates
            while len(split) > 1:
                #Get the Sub-Dictionary Name
                sub = split.pop(0)
                if sub[1:] == '-Channel':
                    sub = int(sub[0])
                
                #Initialize the Sub Dictionary
                if sub not in dictionary:
                    dictionary[sub] = {}
                
                #Set the New Dictionary
                dictionary = dictionary[sub]
            
            #Load the File
            dictionary[split[0].split('.')[0]] = load(filename)
        except:
            pass
    
    #Return the Dictionary
    return gates

#Load the Gates
gates = loadGates()

#Fetch a Gate
def fetch(gate, channel = 2):
    """Fetch a Logic Gate by Name"""
    
    #Get the Gate
    if type(gate) == str:
        if not channel in range(1, 3):
            for i in range(1, 3):
                if gate in gates[i]:
                    gate = gates[i][gate]
                    break
        elif gate in gates[channel]:
            gate = gates[channel][gate]
        else:
            gate = load(gate)
    elif type(gate) == np.ndarray and len(np.shape(gate)) == 2:
        gate = Hashiwokakero(gate)
    elif type(gate) != Hashiwokakero:
        raise TypeError('Unable to cast Type ' + str(type(gate)) + ' to Hashiwokakero')
    
    #Return the Gate
    return gate

#Make a Truth Table
def truthTable(gate, channel:int = 2):
    """Make a Truth Table by Solving a Logic Gate Puzzle
    
    Parameters
    ----------
    gate : str or Hashi
        Either the name of the logic gate, or a Hashi puzzle object for the 
        gate.
    channel : int, optional
        The gate channel to determine which version of the gate to use (and 
        which logic levels are associated with that gate).
        The default is 2.
    
    """
    
    #Get a Fresh Copy of the Gate
    gate = fetch(gate, channel)
    
    #Get the Number of Inputs/Outputs
    M = len(gate.inputs)
    N = len(gate.outputs)
    
    #Print a Header
    print(' | '.join([
        ' '.join(chr(i) for i in range(65, 65 + M)),
        ' '.join(chr(i) for i in range(65 + M, 65 + M + N))]))
    print('|'.join(['--'*M, '--'*N]))
    
    #Compute the Table
    M = 2**M
    for i in range(M, 2*M):
        #Compute the Inputs
        inputs = bin(i)[3:]
        
        #Solve the Puzzle
        outputs = gate.solve(inputs)
        
        #Print the Row in the Truth Table
        print(' | '.join([
            ' '.join(list(inputs)),
            ' '.join(list(outputs))]))

#Time a Gate (Wall Time)
def timeGate(gate, N:int = 10, channel:int = 2, mode:int = 0, Print:bool = True, Return:bool = False):
    """Time a Gate
    
    Parameters
    ----------
    gate : str or Hashiwoakkero
        Either the name of the logic gate, or a Hashiwokakero puzzle object 
        for the gate.
    N : int, optional
        The number of times to test each unique input (the total number of 
        times the gate is solved is equal to N*2**(#inputs)).
        The default is 10.
    channel : int, optional
        The gate channel to determine which version of the gate to use (and 
        which logic levels are associated with that gate).
        The default is 2.
    Print : bool, optional
        Determines if the result of the test should be printed.
        The default is True.
    Return : bool, optional
        When True, all the raw data is returned. The data will also be 
        returned if Print is False (since otherwise there's the computer 
        would just be spinning when this function is run).
        The default is False.
    
    Returns
    -------
    T, uT : float, float
        The average solve time, and the standard deviation in the average 
        solve time.
    
    """
    
    #Get the Gate
    gate = fetch(gate, channel)
    
    #Compute the Number of Inputs
    M = 2**len(gate.inputs)
    
    #Check the Mode
    if mode == 0:
        #Print a Header
        if Print:
            print('Breakdown by Input:')
        
        #Collect Data
        data = np.zeros((M, N))
        averages = np.zeros(M)
        if N > 1:
            stdevs = np.zeros(M)
        for i in range(M):
            #Compute the Inputs
            inputs = bin(M + i)[3:]
            
            #Do N Trials for those Inputs
            for j in range(N):
                #Record the Start Time
                data[i, j] -= time.time()
                
                #Solve the Puzzle
                outputs = gate.solve(inputs)
                
                #Record the Stop Time
                data[i, j] += time.time()
            
            #Compute the Average Solve Time
            averages[i] = np.mean(data[i])
            
            #Compute the Uncertainty on the Mean
            if N > 1:
                stdevs[i] = np.std(data[i], ddof = 1)/np.sqrt(N)
            
            #Print the Result (in ms)
            if Print:
                #Append the Inputs/Outputs
                headder = ''.join((inputs, '|', outputs, ':  '))
                
                #Print the Apropriate Version
                if N == 1:
                    print(str(sfRound(1000*averages[i], 3)).join((headder, ' ms')))
                else:
                    print(' ± '.join(map(str, uRound(1000*averages[i], 1000*stdevs[i]))).join((headder, ' ms')))
    elif mode == 1:
        #Defind a Random Input Generator
        lo = M
        hi = 2*M
        def random():
            return bin(np.random.randint(lo, hi))[3:]
        
        #Time the Gate
        def test():
            #Mark the Starting time
            t = -time.time()
            
            #Solve the Gate
            gate.solve(random())
            
            #Mark the Ending Time
            t += time.time()
            
            #Return the Runtime
            return t
        
        #Generate the Data
        data = np.array([test() for i in range(N)])
    
    #Print the Overall Average
    if Print:
        #Compute the Averare/Uncertainty on the Mean (in ms)
        t = 1000*np.mean(data)
        ut = 1000*np.std(data, ddof = 1)/np.sqrt(N)
        
        #Print the Averate Solve Time
        print('\nAverage Solve Time: %s ± %s ms' % tuple(map(str, uRound(t, ut))))
    
    #Return the Data
    if Return or not Print:
        return data

#Analyze the Runtime of the Solver
def analyze(N:int = 100, Return:bool = True):
    """Time Every Gate in the 'gates' Dictionary, and Plot the Results"""
    #Find/Yield Every Gates
    def Gates(d = gates):
        #Explore the Current Dictionary
        for k, v in d.items():
            #Check the Object
            if type(v) == Hashiwokakero:
                #Yield the Puzzle
                yield k, v
            elif type(v) == dict:
                #Recursively Explore the Sub-Dictionary
                yield from Gates(v)
    
    #Collect Data
    vertices = []
    edges = []
    bridges = []
    ubridges = []
    times = []
    utimes = []
    for name, gate in Gates():
        #Print a Status Update
        print(': '.join(['%s-Channel' % ['0','1','2b','2a'][gate.type], name, '']), end = '')
        
        #Time the Gate
        data = timeGate(gate, N, None, 1, False, True)
        
        #Save the Data
        vertices.append(len(gate.vertices))
        edges.append(len(gate.edges))
        bridges.append((gate.SUM + (len(gate.inputs) + len(gate.outputs))*sum(gate.logic)/2)/2)
        ubridges.append((len(gate.inputs) + len(gate.outputs))*gate.type/2)
        times.append(np.mean(data))
        utimes.append(np.std(data, ddof = 1))
        
        #Finish Printing the Status Update
        try:
            print('%s ± %s ms' % tuple(map(str, uRound(1000*times[-1], 1000*utimes[-1]))))
        except:
            print('%s ms' % str(sfRound(1000*times[-1], 3)))
    
    #Make a Single Data Structure
    data = np.vstack((vertices, edges, bridges, ubridges, times, utimes))
    
    #Plot the Data
    plt.figure(figsize = (4, 8))
    gs = plt.GridSpec(3, 1)
    for i, (xlabel, xerror) in enumerate(zip(['Vertices','Edges','Bridges'],[None, None, data[3]])):
        plt.subplot(gs[i])
        plt.errorbar(data[i], 1000*data[4], 1000*data[5], xerror, '+')
        plt.xlabel(xlabel)
        plt.ylabel('Run Time (ms)')
    plt.show()
    
    #Return the Data
    if Return:
        return data

#Add 2 Numbers Using the X-Bit Adder (this is just for kicks lol)
def add(x:int, y:int, bits:int = 8) -> int:
    """
    Add 2 Numbers Using a Hashiwokakero Puzzle
    
    Parameters
    ----------
    x : int
        The first number you want to add.
    y : int
        THe second number you want to add.
    bits : int, optional
        Which puzzle you want to use to add x and y. Valid 
        options are 2, 4, and 8. Note that your inputs should 
        both be valid integers for the number of bits you've 
        specified (i.e. 2-bit numbers are 0-3, 4-bit numbers 
        are 0-15, and 8-bit numbers are 0-255 inclusive).
        The default is 8.
    
    Returns
    -------
    int
        x + y (computed using the specified puzzle)
    
    """
    
    #Check the Number of Bits
    if bits in (2, 4, 8):
        #Check the Inputs
        if x >= (1 << bits) or x < 0:
            #X is Out of Range
            raise ValueError('x must be an integer on the interval [0, 2**bits). Current value: %d' % x)
        elif y >= (1 << bits) or y < 0:
            #Y is Out of Range
            raise ValueError('y must be an integer on the interval [0, 2**bits). Current value: %s' % y)
        else:
            #Use the Specified Adder
            return int(gates[2][str(bits) + 'b adder'].solve(bin((1 << 2*bits)|(x << bits)^y)[3:]), 2)
    else:
        #Raise an Exception
        raise ValueError('Unsupported number of bits. Supported numbers of bits are 2, 4, and 8.')

#Subtract 2 Numbers Using the 8-Bit Adder (this is also just for kicks)
def subtract(x:int, y:int) -> int:
    """Compute (x - y) Using the 8-Bit 2's Compliment and the 8-Bit Adder"""
    
    #Compute the 2's Compliement of y
    y = gates[2]['8b compliment'].solve(bin(256 + y)[3:])
    
    #Compute x - y
    z = gates[2]['8b adder'].solve(bin(256 + x)[3:] + y)
    
    #Check the Sign
    if z[0] == '1':
        #Return the Number as is
        return int(z[1:], 2)
    else:
        #Return the 2's Compliment of the Current Number (Since it's Negative)
        return -int(gates[2]['8b compliment'].solve(z[1:]), 2)
