from peewee import (FloatField, ForeignKeyField, IntegerField, ManyToManyField,
                    Model, PrimaryKeyField, SqliteDatabase, TextField)


class Genome(Model):
    genome_id = PrimaryKeyField()
    genome_code = TextField()
    genome_name = TextField()

    class Meta:
        table_name = 'genomes_table'

class Contig(Model):
    contig_id = PrimaryKeyField()
    contig_code = TextField()
    genome = ForeignKeyField(Genome, backref='contigs')

    class Meta:
        table_name = 'contigs_table'


class NodeOG(Model):
    node_id = PrimaryKeyField()
    node_name = TextField()
    contig = ForeignKeyField(Contig, backref='nodes')
    description = TextField()
    start_coord = IntegerField()
    end_coord = IntegerField()
    
    class Meta:
        table_name = 'nodes_table'

class NodeKey(Model):
    node_id = PrimaryKeyField()
    node_name = TextField()

class Complexity(Model):
    node = ForeignKeyField(NodeOG, backref='complexity')
    contig = ForeignKeyField(Contig, backref='complexity')
    window_complexity = FloatField()
    prob_window_complexity = FloatField()
    io_complexity = FloatField()
    prob_io_complexity = FloatField()
    window = IntegerField()

    class Meta:
        table_name = 'complexity_table'

class Edge(Model):
    source = ForeignKeyField(NodeKey, backref='edges_sources')
    target = ForeignKeyField(NodeKey, backref='edges_targets')
    frequency = IntegerField()
    genomes = TextField()

    class Meta:
        table_name = 'edges_table'



