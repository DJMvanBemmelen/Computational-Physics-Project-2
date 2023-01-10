import random as rd

class Polymer:
    """
    The Polymer class represents single polymers. It's properties track the positions of each subunit and the weights the polymer had while growing to it's final length.
    """ 
    
    
    def __init__(self):
        self.positions = [(0, 0)]
        self.weight = []
        
    def grow_free(self):
        """
        Updates weights and grows the polymer by one subunit in a random direction, crossing and backfolding allowed
        """
        (end_x, end_y) = self.positions[-1]
        
        north = (end_x, end_y + 1)
        east = (end_x + 1, end_y)
        south = (end_x, end_y - 1)
        west = (end_x - 1, end_y)
        
        directions = [north, east, south, west]

        # update weight
        if len(self.weight) == 0:
            self.weight.append(len(directions))
        else:
            self.weight.append(self.weight[-1] * len(directions))

        # choose randomly from directions
        new_pos = rd.choice(directions)
        self.positions.append(new_pos)

    
    def grow(self):
        """
        Checks available positions, updates the weight and grows the polymer by one subunit in a random possible direction.
        """
        # check available spaces
        (end_x, end_y) = self.positions[-1]
        
        north = (end_x, end_y + 1)
        east = (end_x + 1, end_y)
        south = (end_x, end_y - 1)
        west = (end_x - 1, end_y)
        
        directions = [north, east, south, west]
        possible = []
        
        for position in directions:
            if position in self.positions:
                continue
                
            possible.append(position)
                
        # update weight
        if len(self.weight) == 0:
            self.weight.append(len(possible))
        else:
            self.weight.append(self.weight[-1] * len(possible))
            
        # check if growing is possible
        if len(possible) == 0:
            end_pos = self.positions[-1]
            self.positions.append(end_pos)
            return
        
        # grow the polymer
        new_pos = rd.choice(possible)
        self.positions.append(new_pos)

    
    def prune(self):
        """
        Prune step of the PERM algorithm that either stops the growth of the polymer or multiplyes the weight by 2
        """
        choice = rd.choice([1,2])

        if choice == 1:
            self.weight[-1] = 0
        else:
            self.weight[-1] = self.weight[-1]*2
    

    def enrich_weight(self):
        """
        First part of the Enrichting step of the PERM algortihm that cuts the weight in half
        """
        self.weight[-1] = self.weight[-1]*0.5

    
    def enrich_copy(self, positions, weight):
        """
        Second part of the Enrichting step of the PERM algortihm that copies the polymer 
        """
        for position in positions:
            self.positions.append(position)

        for single_weight in weight:
            self.weight.append(single_weight)
        # self.positions = positions
        # self.weight = weight

    
    def build_single(self, L_max):
        """
        Builds an entire polymer of length L, if no further growth is possible, it breaks off the process.

        Parameters:
        -----------
        L_max : int
            The desired total length of the polymer.
        """
        for i in range(L_max):
            self.grow()
    
    
    def end_to_end(self, L):
        """
        Calculates the square of the end-to-end distance of the polymer.
        
        Parameters:
        -----------
        L: int
            The desired nr of points for the calculation.
        """
        (start_x, start_y) = self.positions[0]
        (end_x, end_y) = self.positions[L]
        
        dist = (end_x - start_x)**2 + (end_y - start_y)**2
        
        return dist
    
    
    def r_gyration(self, nr_points):
        """
        Calculates the square of the radius of gyration of the polymer.
        
        Parameters:
        -----------
        nr_points: int
            The desired nr of points for the calculation.
        """
        # obtain subset for calculation
        x_list = []
        y_list = []
        
        # @@@ Maybe this can be faster ?? @@@
        for i in range(nr_points):
            x_list.append(self.positions[i][0])
            y_list.append(self.positions[i][1])
            
        # calculate c.o.m. position    
        cm_x = 1/nr_points * sum(x_list)
        cm_y = 1/nr_points * sum(y_list)
        
        # calculate radius of gyration
        r_gyration = 0
        
        for i in range(nr_points):
            dist = (x_list[i] - cm_x)**2 + (y_list[i] - cm_y)**2
            r_gyration += dist
            
        r_gyration = 1/nr_points * r_gyration
        
        return r_gyration
        