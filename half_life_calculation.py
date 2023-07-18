#!/usr/bin/env python
# coding: utf-8

# In[53]:


#importing required modules
import pandas as pd
import numpy as np
import warnings
warnings.filterwarnings('ignore')
from sklearn.linear_model import LinearRegression


# In[54]:


#importing Library
import os

#printing the current directory
datapath = os.getcwd()


# In[55]:


##attaching the file to the datapath
file_path = datapath + "\\DecayTimecourse.txt"
file_path


# In[56]:


#Reading the DecayTimecourse text file
raw_data = pd.read_csv(file_path, sep = "\t", skiprows = 1)


# In[57]:


#renaming the first column as "gene_name"
raw_data.rename(columns={ raw_data.columns[0]: "gene_name" }, inplace = True)


# In[58]:


#Total number of unique genes in the dataframe
raw_data.gene_name.nunique()


# In[90]:


#defining a function to fit the data and calculate the half life
def get_half_life(gene_list,df1):
    
    #creating an empty list
    gene_slope_list = []

    for gene in gene_list: #writing a for loop to calculate half life for each gene and concate the resulting dataframe
        df = df1[df1.gene_name == gene] #subsetting for each gene 
        df = df.reset_index(drop = True) 
        time_steps = list(df.columns) #making a list of all columns
        time_steps = time_steps[1:] #list of given time intervals
        gene_exp = df.values.tolist()[0] #converting the gene expression values to a list
        gene_exp = gene_exp[1:] #excluding the first column
        gene_exp_df = pd.DataFrame({'X': time_steps,'y': gene_exp}) #creating a data frame with X being the time intervals and y = gene expression levels
        gene_exp_df = gene_exp_df.dropna(axis = 0) #dropping the time intervals with null values
        gene_exp_df['y'] = gene_exp_df['y'].replace(0, 0.000001) #If the data has 0, it is replaced by 0.000001 as ln(0) is infinity
        if len(gene_exp_df) > 0:
        
            min_val = min(gene_exp_df['y'])
            if min_val <= 0:
                log_contant = 1 - min_val
                gene_exp_df['y'] = gene_exp_df['y'] + log_contant

            gene_exp_df['y'] = np.log(gene_exp_df.y)
            X = gene_exp_df['X'] #defining X
            y = gene_exp_df['y'] #defining y
            X = gene_exp_df.iloc[:, 0].values.reshape(-1, 1)  #values are converted into a numpy array
            y = gene_exp_df.iloc[:, 1].values.reshape(-1, 1) #-1 means that calculate the dimension of rows, but have 1 column

            reg = LinearRegression().fit(X,y) #fitting the values to linear regression model
            reg.intercept_ = reg.intercept_ + 0.00000000001 #to avoid diving by zero

            half_life = np.log(2)/(reg.intercept_) #calculating the t1/2 using formula ln(2)/slope(k)

            save_gene = pd.DataFrame(columns = ['gene','half_life']) #creating a data frame with 2 columns
            save_gene["half_life"] = half_life #one of the column has half_life value
            save_gene["gene"] = str(gene) #another column has the gene name
            gene_slope_list.append(save_gene) 
    t_half_gene = pd.concat(gene_slope_list)
    
    return t_half_gene


# In[91]:


#defining the first data frame
df1 = raw_data.iloc[:,0:10]


# In[92]:


#listing of all the genes from data frame 1
gene_list1 = df1.gene_name.unique()


# In[93]:


#calling the defined function to calculate half life
t_half1_gene = get_half_life(gene_list1, df1)


# In[94]:


#defining the second data frame for df1
df2 = pd.concat([raw_data['gene_name'],raw_data.iloc[:,10:19]], axis = 1)


# In[95]:


#listing of all the genes from data frame 2
gene_list2 = df2.gene_name.unique()


# In[96]:


#calling the defined function to calculate half life
t_half2_gene = get_half_life(gene_list2, df2)


# In[97]:


#defining the third data frame for df2
df3 = pd.concat([raw_data['gene_name'],raw_data.iloc[:,19:28]], axis = 1)


# In[98]:


#listing of all the genes from data frame 3
gene_list3 = df3.gene_name.unique()


# In[99]:


#calling the defined function to calculate half life for df3
t_half3_gene = get_half_life(gene_list3, df3)


# In[100]:


#merging the data frames of first and second time course half life
mer_col = pd.merge(t_half1_gene, t_half2_gene, on="gene", how="outer")


# In[101]:


#merging the above dataframe with the third time course half life 
final = pd.merge(mer_col, t_half3_gene, on="gene", how="outer")


# In[102]:


final.head(10)


# In[103]:


#calculating the average of obtained half lives
final['avg_half_life'] = final.mean(axis=1)


# In[104]:


final


# In[105]:


#sorting the data frame based average half life in ascending order
sorted_df = final.sort_values(by=['avg_half_life'], ascending=True)
sorted_df


# In[111]:


ten_percent_count = int(len(sorted_df)*.1)


# In[112]:


#the genes with very low half lives (bottom 10%)
bottom_10_percent = sorted_df.head(ten_percent_count)


# In[113]:


bottom_10_percent


# In[114]:


#the genes with very high half lives (top 10%)
top_10_percent = sorted_df.tail(ten_percent_count)


# In[115]:


top_10_percent

