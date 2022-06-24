## Species specific rasters


## Set working environment ####
rm(list = ls())
gc()

x = c("DBI", "RPostgreSQL", "data.table", "rpostgis", "sp")
lapply(x, require, character.only = TRUE)
rm(x)

## Connect to server
source("~/gsdms_r_vol/tempdata/workdir/gbifprocessing/scripts/connect_to_server.R")

## Check db has PostGIS installed
pgPostGIS(con)
pgListGeom(con, geog = TRUE)



x <- dbGetQuery(con, "
            SELECT ST_AsText( t1.wkb_geometry )
            FROM public.ecoregions_dinerstein_2017 t1, public.gbif_reptilia t2
            WHERE t2.species = 'Trachemys scripta'
            AND ST_Intersects(t1.wkb_geometry, t2.points_geom);")
saveRDS(y, file = file.path("~/gsdms_r_vol/tempdata/workdir/gsdms/temp/stastext.rds"))

y <- dbGetQuery(con, "
            SELECT ST_AsEWKT( t1.wkb_geometry )
            FROM public.ecoregions_dinerstein_2017 t1, public.gbif_reptilia t2
            WHERE t2.species = 'Anguis fragilis'
            AND ST_Intersects(t1.wkb_geometry, t2.points_geom);")



## DO it other way around: poly > get all points > label points with poly id

SELECT t2.* 
FROM public.ecoregions_dinerstein_2017 t1, public.gbif_reptilia t2
WHERE geom && (SELECT geom FROM boundaries WHERE t1.eco_name = "Adelie Land tundra");










-- SELECT t1.wkb_geometry
-- FROM public.ecoregions_dinerstein_2017 t1, public.gbif_reptilia t2
-- WHERE t2.species = 'Trachemys scripta'
-- AND ST_Contains(t1.wkb_geometry, t2.points_geom);



--    AND 
--   ST_Intersects(ST_Buffer(t2.pt,t2.radius),t1.area)


--   SELECT (wkb_geometry) FROM
-- public.ecoregions_dinerstein_2017 as A AND public.gbif_reptilia AS B
-- WHERE A.wkb_geometry


-- SELECT a.*
-- FROM public.ecoregions_dinerstein_2017 a, public.gbif_reptilia 
-- WHERE EXISTS (SELECT b.geom 
-- FROM b 
-- WHERE b.species == 'Trachemys scripta'
-- AND ST_Intersects(a.geom, b.geom));


-- SELECT ST_INTERSECTION(a.geom, b.geom), 'fair'
-- FROM public.ecoregions_dinerstein_2017 a, public.gbif_reptilia b
-- WHERE a.ID < b.ID
-- AND ST_INTERSECTS(a.geom, b.geom);


-- SELECT id, geom FROM polygon
-- WHERE ST_Intersects(geom, 'SRID=4326;POINT(1 2)')
-- ORDER BY id
-- LIMIT 1000 OFFSET 0;
