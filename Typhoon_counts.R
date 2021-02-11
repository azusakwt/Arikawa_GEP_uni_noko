# 2017
# July:
#   TS3, TS5,
# Aug:
#   None
# 

# 2018
# July:
#   TS7, TS12
# Aug:
#   TS12, TS15, TS19, TS21
# 

# 2019
# July:
#   TS5
# Aug:
#   TS8, TS10
# 

# 2020
# July: 
#   None
# Aug: 
#   TS5, TS9

tibble(year = 2017:2020,
       Jul = c(2,2,1,0),
       Aug = c(0,4,2,2)) %>% 
  mutate(total = Jul + Aug)
