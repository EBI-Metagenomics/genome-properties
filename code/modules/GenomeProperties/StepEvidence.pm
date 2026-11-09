package GenomeProperties::StepEvidence;

sub new {
  my $class = shift;
  my $hashRef = shift;  
  my $self = {};
   #Need to check that this does not go out of scope! 
  $self=$hashRef;

  bless( $self, $class);
  
  return $self;

}

sub type {
  my ($self) = @_;

  return $self->{type};

}

sub sufficient {
  my ($self) = @_;
  return (defined($self->{type}) and $self->{type} eq 'sufficient') ? 1 : 0;
}

#Formerly HMM-CLUST: only counts when the gene is near the rest of the property.
sub clustered {
  my ($self) = @_;
  return (defined($self->{type}) and $self->{type} eq 'clustered') ? 1 : 0;
}

sub accession {
  my ($self) = @_;  
  return $self->{accession};
}


sub step_ev_id {
  my ($self) = @_;  
  return $self->{step_ev_id};
}


sub step {
  my ($self) = @_;  
  return $self->{step};
}


sub get_go {
  my ($self) = @_;  
  return $self->{get_go};
}

sub interpro {
  my ($self) = @_;  
  return $self->{interpro};
}

sub gp {
  my ($self) = @_;  
  return $self->{genome_property};
}

1;
