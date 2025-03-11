# -*- coding: utf-8 -*-
"""
Created on Mon Mar 10 14:27:19 2025

@author: fabio
"""

class Reorder:
    def __init__(self, list1, list2):
        self.list1 = list1
        self.list2 = list2
        self.order = ["CH4",	  "C2H6",	  "C3H8",	  "iC4H10",  "nC4H10",  "iC5H12",  "nC5H12",  "nC6H14",  "nC7H16",  "nC8H18",	
                       "nC9H20",  "nC10H22", "nC11H24", "nC12H26", "nC14H30", "N2", "H2O", "CO2", "C15+"]
    
    def reorder(self):
        combined = dict(zip(self.list1, self.list2))
        ordered_pairs = [(key, combined[key]) for key in self.order if key in combined]
        
       # if not ordered_pairs:
       #     raise ValueError("Nenhuma correspondência encontrada entre os nomes fornecidos e a ordem desejada.")
        
        new_list1, new_list2 = zip(*ordered_pairs)
        return list(new_list1), list(new_list2)

               
               
    
    
    