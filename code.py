import numpy as np
import matplotlib.pyplot as plt


# NO other imports are allowed

class Shape:
    '''
    DO NOT MODIFY THIS CLASS

    DO NOT ADD ANY NEW METHODS TO THIS CLASS
    '''

    def __init__(self):
        self.T_s = None
        self.T_r = None
        self.T_t = None

    def translate(self, dx, dy):
        '''
        Polygon and Circle class should use this function to calculate the translation
        '''
        self.T_t = np.array([[1, 0, dx], [0, 1, dy], [0, 0, 1]])

    def scale(self, sx, sy):
        '''
        Polygon and Circle class should use this function to calculate the scaling
        '''
        self.T_s = np.array([[sx, 0, 0], [0, sy, 0], [0, 0, 1]])

    def rotate(self, deg):
        '''
        Polygon and Circle class should use this function to calculate the rotation
        '''
        rad = deg * (np.pi / 180)
        self.T_r = np.array([[np.cos(rad), np.sin(rad), 0], [-np.sin(rad), np.cos(rad), 0], [0, 0, 1]])

    def plot(self, x_dim, y_dim):
        '''
        Polygon and Circle class should use this function while plotting
        x_dim and y_dim should be such that both the figures are visible inside the plot
        '''
        x_dim, y_dim = 1.2 * x_dim, 1.2 * y_dim
        plt.plot((-x_dim, x_dim), [0, 0], 'k-')
        plt.plot([0, 0], (-y_dim, y_dim), 'k-')
        plt.xlim(-x_dim, x_dim)
        plt.ylim(-y_dim, y_dim)
        plt.grid()
        plt.show()


class Polygon(Shape):
    '''
    Object of class Polygon should be created when shape type is 'polygon'
    '''

    def __init__(self, A):
        '''
        Initializations here
        '''
        self.A = A                                                                        # initialising the standard matrix
        Shape.__init__(self)                                                              # calling init from parent class
        self.original = 0                                                                 #creating a copy of standard matrix to help with plotting

    def translate(self, dx, dy):
        '''
        Function to translate the polygon

        This function takes 2 arguments: dx and dy

        This function returns the final coordinates
        '''
        Shape.translate(self, dx, dy)                                                       # calling translate from parent class
        self.original = (self.A).copy()                                                     #Creating a copy to help with plotting
    
        
        
        for i in range(len(self.A)):
            self.A[i] = np.dot(self.T_t, [self.A[i][0], self.A[i][1], 1])                   # performing translation transformation to get updated coordinates

        x2, y2 = [], []
        for i in range(len(self.A)):
            x2.append(self.A[i][0])
            y2.append(self.A[i][1])

        x2 = (np.around(np.array(x2), 2))                                                   # Rounding off x and y coordinates to 2 decimal places
        y2 = (np.around(np.array(y2), 2))

        return x2, y2

    def scale(self, sx, sy):                                                            
        '''
        Function to scale the polygon

        This function takes 2 arguments: sx and sx

        This function returns the final coordinates
        '''
        Shape.scale(self, sx, sy)                                                           # Calling scale from parent class
        self.original = (self.A).copy()
        

        x1, y1 = [], []

        for i in range(len(self.A)):
            x1.append(self.A[i][0])
            y1.append(self.A[i][1])

        mid_x = sum(x1) / (len(x1))                                                         # calculating the mid points for scalings
        mid_y = sum(y1) / (len(y1))                                                         
        for i in range(len(self.A)):
            self.A[i] = np.dot(self.T_s, [self.A[i][0] - mid_x, self.A[i][1] - mid_y, 1])   #perforing rounding trasnformation to get updated coordinates

        x2, y2 = [], []
        for i in range(len(self.A)):
            x2.append(self.A[i][0] + mid_x)                                                 # updating the final coordinates by adding the centre point
            y2.append(self.A[i][1] + mid_y)
        x2 = (np.around(np.array(x2), 2))
        y2 = (np.around(np.array(y2), 2))

        for i in range(len(self.A)):
            self.A[i][0] += mid_x                                                           # Updating the A matrix
            self.A[i][1] += mid_y

        return x2, y2

    def rotate(self, deg, rx=0, ry=0):
        '''
        Function to rotate the polygon

        This function takes 3 arguments: deg, rx(optional), ry(optional). Default rx and ry = 0. (Rotate about origin)

        This function returns the final coordinates
        '''
        Shape.rotate(self, deg)                                                             # Calling rotation from parent class
        self.original = (self.A).copy()

       
        for i in range(len(self.A)):
            self.A[i] = np.dot(self.T_r, [self.A[i][0] - rx, self.A[i][1] - ry, 1])         #perforing rotation trasnformation to get updated coordinates

        x2, y2 = [], []
        for i in range(len(self.A)):
            x2.append(self.A[i][0] + rx)                                                    # updating the final coordinates by adding the points to be rotated about
            y2.append(self.A[i][1] + ry)
       

        x2 = (np.around(np.array(x2), 2))
        y2 = (np.around(np.array(y2), 2))
        return x2, y2

    def plot(self):
        '''
        Function to plot the polygon

        This function should plot both the initial and the transformed polygon

        This function should use the parent's class plot method as well

        This function does not take any input

        This function does not return anything
        '''

        x1, y1 = [], []

        for i in range(len(self.original)):                                         #loop to plot the initial polygon
            x1.append(self.original[i][0])
            y1.append(self.original[i][1])
        x1.append(self.original[0][0])
        y1.append(self.original[0][1])

        x2, y2 = [], []
        for i in range(len(self.A)):                                                #loop to plot the final polygon
            x2.append(self.A[i][0])
            y2.append(self.A[i][1])
        x2.append(self.A[0][0])
        y2.append(self.A[0][1])
        
        plt.plot(x1, y1, 'k', ls=':')                                               #plotting the initial matrix with dashed line
        plt.plot(x2, y2, 'k')                                                       #plotting the final matrix 
        
            
        x2.extend(x1)
        y2.extend(y1)
        x_dim = max(map(abs, x2))                                                   #setting the limits as the max of x and y
        y_dim = max(map(abs, y2))
        Shape.plot(self, x_dim, y_dim)                                              #calling plot function from parent class


class Circle(Shape):
    '''
    Object of class Circle should be created when shape type is 'circle'
    '''

    def __init__(self, x=0, y=0, radius=5):
        '''
        Initializations here
        '''
        self.x = x
        self.y = y                                                               #Initialising the x and y coordinates
        self.x_og=0.0                                                           #initialising temp variables to store original coordinates
        self.y_og=0.0
        self.radius = radius                                                    #Initialising the radius of teh circle
        self.og_radius = 0.0                                                    #initialising the temp variable to store original radius
        
        Shape.__init__(self)

    def translate(self, dx, dy):
        '''
        Function to translate the circle

        This function takes 2 arguments: dx and dy (dy is optional).

        This function returns the final coordinates and the radius
        '''
        Shape.translate(self, dx, dy)                                           #calling translate function from parent class
        self.x_og = self.x
        self.y_og = self.y
        self.og_radius=self.radius

        new_coordinates = np.dot(self.T_t, [[self.x], [self.y], [1]])           #perforing trasnlation trasnformation to get updated coordinates


        self.x = new_coordinates[0]
        self.y = new_coordinates[1]

        self.x = round(float(self.x), 2)                                        #rounding the x and y coordinated to 2 deciaml places
        self.y = round(float(self.y), 2)
        self.radius = round(float(self.radius), 2)

        return self.x, self.y, self.radius                                      #returning updated x y coordinates and radius

    def scale(self, sx):
        '''
        Function to scale the circle

        This function takes 1 argument: sx

        This function returns the final coordinates and the radius

        '''
        self.x_og = self.x
        self.y_og = self.y
        self.og_radius=self.radius
        Shape.scale(self, sx, 0)                                                #calling scale from the parent class
       
        new_coordinates = np.dot(self.T_s, [[self.radius], [self.radius], [1]])  #tranforming the matrix to scale
        self.radius = new_coordinates[0][0]                                     #updating the radius

        

        self.x = round(float(self.x), 2)                                             #rounding the x and y coordinated to 2 deciaml places
        self.y = round(float(self.y), 2)
        self.radius = round(self.radius, 2)
        return self.x, self.y, self.radius                                       #returning the updated x y and radius coordinates

    def rotate(self, deg, rx=0, ry=0):
        '''
        Function to rotate the circle

        This function takes 3 arguments: deg, rx(optional), ry(optional). Default rx and ry = 0. (Rotate about origin)

        This function returns the final coordinates and the radius
        '''
        Shape.rotate(self, deg)                                                 #calling the rotate function from the parent class
        self.x_og = self.x                                                      
        self.y_og = self.y
        self.og_radius=self.radius
        
        new_coordinates = np.dot(self.T_r, [(self.x - rx), (self.y - ry), 0])   #transforming the matrix to rotate
        

        self.x = new_coordinates[0] + rx
        self.y = new_coordinates[1] + ry
        plt.axis('scaled')
        self.x = round(self.x, 2)
        self.y = round(self.y, 2)
        self.radius = round(self.radius, 2)
        return self.x, self.y, self.radius                                      #returning the updated x y and radius

    def plot(self):
        '''
        Function to plot the circle

        This function should plot both the initial and the transformed circle

        This function should use the parent's class plot method as well

        This function does not take any input

        This function does not return anything

        '''

        grid = plt.gca()

        circle = plt.Circle((self.x_og, self.y_og), self.og_radius, fill=False, ls=':')    #plotting the initial circle
        new_circle = plt.Circle((self.x, self.y), self.radius, fill=False)                  #plotting the final matrix
        
        
        grid.add_patch(circle)                                                              #figure is added to the plot
        grid.add_patch(new_circle)
       
            

        a = abs(self.x) + self.radius
        b = abs(self.y) + self.radius
        x_dim = max(a, abs(self.x_og)+self.radius)                                      #setting the dimension as the maximum possible value from the given coordinates
        y_dim = max(b, abs(self.y_og)+self.radius)
        Shape.plot(self, x_dim, y_dim)                                              #calling plot from parent class
        


if __name__ == "__main__":
    '''
    Add menu here as mentioned in the sample output section of the assignment document.
    '''

    verbose = int(input('verbose? 1 to plot, 0 otherwise: '))
    n = int(input('Enter the number of test cases :'))
    for i in range(n):
        choice = int(input('which shape would you like to work upon: (1 for Circle , 0 for polygon) :'))

        x_list, y_list = [], []
        if choice == 0 and verbose == 0:
            n = int(input('Enter the number of sides: '))
            for i in range(n):
                x, y = input('Enter x y coordinates: ').split()
                x_list.append(float(x))
                y_list.append(float(y))
            a = np.array([x_list, y_list, [1] * n])
            a = a.transpose()
            obj = Polygon(a)

            q = int(input('Enter number of queries : '))
            print('1) R deg (rx) (ry) \n 2) T dx (dy) \n 3) S sx (sy) \n 4) P')
            for _ in range(q):
                query = input('Enter query as per the format: ').split()
                if query[0] == 'T':
                    dx = query[1]
                    if len(query) == 3:
                        dy = query[2]
                    else:
                        dy = dx

                    x, y = obj.translate(float(dx), float(dy))
                    x1 = x.tolist()
                    y1 = y.tolist()
                    for i in y1:
                        x1.append(i)
                    print(' '.join(list(map(str, x1))))

                if query[0] == 'S':
                    sx = query[1]
                    if len(query) == 3:
                        sy = query[2]
                    else:
                        sy = sx
                    x, y = obj.scale(float(sx), float(sy))
                    x1 = x.tolist()
                    y1 = y.tolist()
                    for i in y1:
                        x1.append(i)
                    print(' '.join(list(map(str, x1))))

                if query[0] == 'R':
                    deg = query[1]
                    if len(query) == 2:
                        x, y = obj.rotate(int(deg))
                        x1 = x.tolist()
                        y1 = y.tolist()
                        for i in y1:
                            x1.append(i)
                        print(' '.join(list(map(str, x1))))
                    elif len(query) == 3:
                        rx = query[2]

                        x, y = obj.rotate(int(deg), float(rx))
                        x1 = x.tolist()
                        y1 = y.tolist()
                        for i in y1:
                            x1.append(i)
                        print(' '.join(list(map(str, x1))))
                    elif len(query) == 4:
                        rx = query[2]
                        ry = query[3]

                        x, y = obj.rotate(int(deg), float(rx), float(ry))
                        x1 = x.tolist()
                        y1 = y.tolist()
                        for i in y1:
                            x1.append(i)
                        print(' '.join(list(map(str, x1))))
                elif query[0] == 'P':
                    obj.plot()

        if verbose == 1 and choice == 0:
            n = int(input('Enter the number of sides: '))
            for i in range(n):
                x, y = input('Enter x y coordinates: ').split()
                x_list.append(float(x))
                y_list.append(float(y))
            a = np.array([x_list, y_list, [1] * n])
            a = a.transpose()
            obj = Polygon(a)

            q = int(input('Enter number of queries : '))
            print('1) R deg (rx) (ry) \n 2) T dx (dy) \n 3) S sx (sy) \n 4) P')
            for _ in range(q):
                query = input('Enter query as per the format: ').split()
                if query[0] == 'T':
                    dx = query[1]
                    if len(query) == 3:
                        dy = query[2]
                    else:
                        dy = dx

                    x, y = obj.translate(float(dx), float(dy))
                    x1 = x.tolist()
                    y1 = y.tolist()
                    for i in y1:
                        x1.append(i)
                    print(' '.join(list(map(str, x1))))
                    obj.plot()

                if query[0] == 'S':
                    sx = query[1]
                    if len(query) == 3:
                        sy = query[2]
                    else:
                        sy = sx
                    x, y = obj.scale(float(sx), float(sy))
                    x1 = x.tolist()
                    y1 = y.tolist()
                    for i in y1:
                        x1.append(i)
                    print(' '.join(list(map(str, x1))))
                    obj.plot()
                    
                if query[0] == 'R':
                    deg = query[1]
                    if len(query) == 2:
                        x, y = obj.rotate(int(deg))
                        x1 = x.tolist()
                        y1 = y.tolist()
                        for i in y1:
                            x1.append(i)
                        print(' '.join(list(map(str, x1))))
                        obj.plot()
                    elif len(query) == 3:
                        rx = query[2]

                        x, y = obj.rotate(int(deg), float(rx))
                        x1 = x.tolist()
                        y1 = y.tolist()
                        for i in y1:
                            x1.append(i)
                        print(' '.join(list(map(str, x1))))
                        obj.plot()
                    elif len(query) == 4:
                        rx = query[2]
                        ry = query[3]

                        x, y = obj.rotate(int(deg), float(rx), float(ry))
                        x1 = x.tolist()
                        y1 = y.tolist()
                        for i in y1:
                            x1.append(i)
                        print(' '.join(list(map(str, x1))))
                        obj.plot()
                elif query[0] == 'P':
                    obj.plot()

        if choice == 1 and verbose == 1:
            z = list(map(float,input('Enter centre cordinates and radius :').split()))
            obj = Circle(z[0], z[1], z[2])

            q = int(input('Enter number of queries : '))
            print('1) R deg (rx) (ry) \n 2) T dx (dy) \n 3) S sx (sy) \n 4) P')
            for _ in range(q):
                query = input('Enter query as per the format: ').split()
                if query[0] == 'T':
                    dx = query[1]
                    if len(query) == 2:
                        dy = dx
                    else:
                        dy = query[2]
                    x, y, radius = obj.translate(float(dx), float(dy))
                    print(x, y, radius)
                    obj.plot()
                    
                if query[0]=='S':
                    sx=query[1]
                    x, y, radius = obj.scale(float(sx))
                    print(x,y,radius)
                    obj.plot()
                    
                if query[0]=='R':
                    deg=query[1]
                    rx,ry=0,0
                    if len(query)==3:
                        rx=query[2]
                       
                    elif len==4:
                        rx=query[2]
                        ry=query[3]
                    x,y,radius=obj.rotate(float(deg),float(rx),float(ry))
                    print(x,y,radius)
                    obj.plot()
                        
                if query[0]=='P':
                    obj.plot()



        if choice == 1 and verbose == 0:
            z = list(map(float,input('Enter centre cordinates and radius :').split()))
            obj = Circle(z[0], z[1], z[2])

            q = int(input('Enter number of queries : '))
            print('1) R deg (rx) (ry) \n 2) T dx (dy) \n 3) S sx (sy) \n 4) P')
            for _ in range(q):
                query = input('Enter query as per the format: ').split()
                if query[0] == 'T':
                    dx = query[1]
                    if len(query) == 2:
                        dy = dx
                    else:
                        dy = query[2]
                    x, y, radius = obj.translate(float(dx), float(dy))
                    print(x, y, radius)
                   
                if query[0]=='S':
                    sx=query[1]
                    x, y, radius = obj.scale(float(sx))
                    print(x,y,radius)
                   
                    
                if query[0]=='R':
                    deg=query[1]
                    rx,ry=0,0
                    if len(query)==3:
                        rx=query[2]
                       
                    elif len==4:
                        rx=query[2]
                        ry=query[3]
                    x,y,radius=obj.rotate(float(deg),float(rx),float(ry))
                    print(x,y,radius)
                  
                        
                if query[0]=='P':
                    obj.plot()
                
                    
                    
