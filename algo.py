# This is an algorithm that balances chemical equations (edge cases to be addressed).
# 
# 1. Parse the equation into reactants and products
# 2. Identify all elements
# 3. Build a system of linear equations from the chemical equation
# 4. Find the null space of the homegenous matrix by rref
# 5. Scale properly to take care of fractional solutions

# Start with a parsing algo 

import sympy
from math import gcd

class ChemicalEquationBalancer():
    def __init__(self):
        self.elements = {} # Elements with their qunatities
        self.compounds = [] # All compounds from reactants and products
    
    def parse_equation(self, equation):
        # Remove the spaces
        equation = equation.replace(" ", "").replace("_", "").replace("\longrightarrow", "-->")

        # Split reactants and products
        sides = equation.split("-->")
        if len(sides) != 2:
            raise ValueError("Unable to parse the equation")
        
        # Split compunds
        reactants = sides[0].split("+")
        products = sides[1].split("+")

        self.compounds = reactants + products

        return reactants, products
    
    def _parse_compound(self, compound):
        """Extract elements in a compound and their quantities"""
        elements = {}

        # All chemical elements have at most one lowercase letter and at most two letters
        # Compounds may have parentheses, in which there are subcompounds
        i = 0;
        while i < len(compound):
            # If there are subcompounds
            if compound[i] == '(':
                j = i + 1
                paren_count = 1
                while paren_count > 0 and j < len(compound):
                    if compound[j] == '(':
                        paren_count += 1
                    elif compound[j] == ')':
                        paren_count -= 1
                    j += 1

                subcompound = compound[i + 1: j - 1]
                sub_elements = self._parse_compound(subcompound) # returns a dictionary

                # Check for a subscript
                subscript = ""
                while j < len(compound) and compound[j].isdigit():
                    subscript += compound[j]
                    j += 1
                
                count = int(subscript) if subscript else 1

                for element, quantity in sub_elements.items():
                    if element in elements:
                        elements[element] += quantity * count
                    else:
                        elements[element] = quantity * count
                i = j

            # If there are no subcompounds
            elif compound[i].isupper():
                if i + 1 < len(compound) and compound[i + 1].islower():
                    element = compound[i: i + 2]
                    i += 2
                else:
                    element = compound[i]
                    i += 1

                # Check for a subscript
                subscript = ""
                while i < len(compound) and compound[i].isdigit():
                    subscript += compound[i]
                    i += 1  

                count = int(subscript) if subscript else 1  

                if element in elements:
                    elements[element] += count
                else:
                    elements[element] = count 
            # Skip certain characters (that come with MathLive rendering)
            else:
                i += 1
        return elements  


    def _build_matrix(self, reactants, products):
        # Collect all elements from all compounds
        all_elements = set()

        # Process reactants
        for compound in reactants:
            elements = self._parse_compound(compound)
            for element in elements:
                all_elements.add(element)

        # Process products
        for compound in products:
            elements = self._parse_compound(compound)
            for element in elements:
                all_elements.add(element)

        # Initialize the matrix with zeros
        # Rows = elements; Columns = compounds
        element_list = list(all_elements)
        matrix = []

        # Build by rows
        for element in element_list:
            row = []

            for compound in reactants:
                elements = self._parse_compound(compound)
                count = -elements.get(element, 0) # negative coefficient
                row.append(count)

            for compound in products:
                elements = self._parse_compound(compound)
                count = elements.get(element, 0)
                row.append(count)
            
            matrix.append(row)

        return matrix, element_list
    
    
    def _solve_matrix(self, matrix):
        # import numpy as np
        # from fractions import Fraction

        # A = np.array(matrix, dtype = float)

        # u, s, vh = np.linalg.svd(A)

        # null_space = vh[-1]

        # fractions = [Fraction(val).limit_denominator() for val in null_space]

        # from math import lcm
        # denominators = [f.denominator for f in fractions]
        # lcm_value = 1
        # for d in denominators:
        #     lcm_value = lcm(lcm_value, d)

        # coefficients = [int(f * lcm_value) for f in fractions]
    
        # if min(coefficients) < 0:
        #     coefficients = [-c for c in coefficients]
        
        # # Simplify by GCD
        # from math import gcd
        # gcd_value = coefficients[0]
        # for c in coefficients[1:]:
        #     gcd_value = gcd(gcd_value, c)
        
        # coefficients = [c // gcd_value for c in coefficients]

        # return coefficients

        # Create a sympy Matrix and compute its nullspace
        M = sympy.Matrix(matrix)
        ns = M.nullspace()
        if not ns:
            raise Exception("No solution found for the matrix.")
        sol = ns[0]
    
        # Multiply by the LCM of denominators to get integer coefficients
        denominators = [x.q for x in sol]  # x.q gives the denominator of a Rational
        lcm_denom = 1
        for d in denominators:
            lcm_denom = sympy.ilcm(lcm_denom, d)
        coeffs = [int(x * lcm_denom) for x in sol]
    
        # Ensure all coefficients are positive
        if min(coeffs) < 0:
            coeffs = [-c for c in coeffs]
    
        # Simplify coefficients by dividing by their GCD
        g = coeffs[0]
        for c in coeffs[1:]:
            g = gcd(g, c)
        coeffs = [c // g for c in coeffs]
    
        return coeffs
    

    def get_balanced_equation(self, coefficients, reactants, products):
        """Format the balanced equation with coefficients"""
        # Format reactants
        balanced_reactants = []
        for i, compound in enumerate(reactants):
            coef = coefficients[i]
            if coef == 1:
                balanced_reactants.append(compound)
            else:
                balanced_reactants.append(f"{coef}{compound}")
        
        # Format products
        balanced_products = []
        for i, compound in enumerate(products):
            coef = coefficients[i + len(reactants)]
            if coef == 1:
                balanced_products.append(compound)
            else:
                balanced_products.append(f"{coef}{compound}")
        
        # Combine into final equation
        balanced_equation = " + ".join(balanced_reactants) + r"\to" + " + ".join(balanced_products)
        
        return balanced_equation
    
    
    def balance(self, equation):
        """Main method to balance a chemical equation"""
        # Parse the equation
        reactants, products = self.parse_equation(equation)
        
        # Build the coefficient matrix
        matrix, elements = self._build_matrix(reactants, products)
        
        # Solve for the coefficients
        coefficients = self._solve_matrix(matrix)
        
        # Format the balanced equation
        balanced_equation = self.get_balanced_equation(coefficients, reactants, products)
        
        return balanced_equation
    
# test_equations = [
#     "H2 + O2 --> H2O",
#     "C + O2 --> CO2",
#     "C6H12O6 + O2 --> CO2 + H2O",
#     "Fe + O2 --> Fe2O3",
#     "KMnO4 + HCl --> KCl + MnCl2 + H2O + Cl2"
# ]
# --- Won't Work due to parsing algo change --- 
# if __name__ == "__main__":
#     balancer = ChemicalEquationBalancer()
#     for eq in test_equations:
#         print(f"\nBalancing: {eq}")
#         try:
#             result = balancer.balance(eq)
#             print(f"Result: {result}")
#         except Exception as e:
#             print(f"Error: {e}")