require(RPostgreSQL)

drv <- dbDriver("PostgreSQL")  #Specify a driver for postgreSQL type database
con <- dbConnect(drv, dbname="qaeco_spatial", user="qaeco", password="Qpostgres15", host="boab.qaeco.com", port="5432")  #Connection to database server on Boab

roads <- dbGetQuery(con,"
                    SELECT
                    r.uid as uid, r.length/1000 AS length, ST_X(r.geom) AS x, ST_Y(r.geom) AS y
                    FROM
                    (SELECT
                    uid, ST_Length(geom) AS length, ST_ClosestPoint(geom, ST_Centroid(geom)) AS geom
                    FROM
                    gis_victoria.vic_gda9455_roads_state) AS r
                    ")

roads2 <- dbGetQuery(con,"
                    SELECT
                    uid, ST_Length(geom) AS length, ST_ClosestPoint(geom, ST_Centroid(geom)) AS geom
                    FROM
                    gis_victoria.vic_gda9455_roads_state
                    ")


gbif <- data.table(dbGetQuery(con,"
                    SELECT *
                    FROM gbif.raw_data
                   WHERE family = 'Felidae'
                   LIMIT 10;
                     "))



gbif[countrycode == "", countrycode := "AU"]

dbWriteTable(con, c("gbif", "temp"), value = gbif, row.names=FALSE, overwrite=TRUE)




