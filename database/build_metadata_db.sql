create or replace table sra_metadata as select * from read_parquet('*') order by acc;
alter table sra_metadata add primary key (acc);
