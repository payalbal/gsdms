
## Connect to server

drv <- dbDriver("PostgreSQL") ## specify driver for postgreSQL database
con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15", host="boab.qaeco.com", port="5432")