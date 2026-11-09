package GenomeProperties;	

use strict;
use warnings;
use Data::Printer;
use Carp;
use GenomeProperties::Definition;
use GenomeProperties::Step;
use GenomeProperties::StepEvidence;
use GenomePropertiesIO;
use Text::Wrap;
use LWP::Simple;
use JSON;
use IO::Uncompress::Gunzip qw($GunzipError);
use List::Util qw(min max);

sub new {
  my $class = shift;
  
  my $self = {};
 
  #These defaults can/should be over-ridden. 
  $self->{gpdir}          = "/nfs/production/interpro/genome_properties/data/genome_properties_flatfiles/";
  $self->{hmm_dir}        = "/nfs/production/interpro/genome_properties/data/pfam_and_tigrfam/";
  $self->{project_name}   = "DEFAULT";
  $self->{debug}          = 0;
  $self->{read_sig}       = 0;
  #Clustered evidence (HMM-CLUST in the JCVI release) only counts when the
  #gene is within this many base pairs of a gene found by other evidence.
  $self->{cluster_distance} = 3000;
  #strict: clustered evidence needs gene positions (-gff), otherwise it is ignored.
  #permissive: clustered evidence is treated like any other evidence.
  #If not set, strict is used when gene positions are given, otherwise permissive.

  #These should be constants;
  $self->{_skip}          = { CATEGORY => 1, ROOT => 1, SUMMARY => 1 };
  
  bless( $self, $class);  
  return $self;
}

sub annotated {
  my $self = shift;
  my $state = shift;

  if(defined($state)){
    $self->{annotated} = 1;
  }
  return ($self->{annotated});
}


sub write_results{
  my ( $self ) = @_;
  #Only the table file has a header line
  if($self->tableFH){
    my $fh = $self->tableFH;
    print $fh  "acc\tdescription\tproperty_value\tstep_number\tstep_name\tstep_value\trequired?\tevidence_type\tevidence_name\tbest_HMM_hit\tHMM_hit_count\n";
  }
  
  foreach my $acc (sort{$a cmp $b}keys(% { $self->get_defs })){
    #next unless($acc eq 'GenProp0643' or $acc eq 'GenProp0001' or $acc eq 'GenProp0639');
    my $prop = $self->get_def($acc);
    next if(defined( $self->{_skip}->{ $prop->type } ));
    $self->print_summary($prop) if($self->summaryFH and fileno($self->summaryFH));
    $self->print_long($prop) if($self->longFH);
    $self->print_table($prop) if($self->tableFH);
    $self->print_compressed_table($prop) if($self->compTableFH);
    $self->build_json($prop) if($self->jsonFH);
    $self->print_matches($prop) if($self->matchFH);
  }

  $self->print_json if($self->jsonFH);
  return(1);
}


sub build_json {
  my ( $self, $prop ) = @_;
  my %propj = ( property => $prop->accession,
                name     => $prop->{name},
                values   => { $self->{name} => $prop->result,
                              TOTAL       => { YES     => $prop->result eq 'YES'     ? 1 : 0,
                                               PARTIAL => $prop->result eq 'PARTIAL' ? 1 : 0,
                                               NO      => $prop->result eq 'NO'      ? 1 : 0      }});


  foreach my $step (sort { $a->order cmp $b->order} @{ $prop->get_steps }){
    next if($step->skip);

    push(@{ $propj{steps} }, { "step"       => $step->order,
                               "step_name"  => $step->step_name,
                               "required"   => defined($step->required) ? 1 : 0,
                               "values"     => { $self->{name} => $step->found } });
  }

  $self->{_json}->{$prop->accession}= \%propj;
  return;
}

sub print_json {
  my ($self) = shift;

  my $jfh = $self->jsonFH;
  print $jfh to_json($self->{_json}, { ascii => 1, pretty => 1 });
  #print $self->jsonFH $j;
  return;
}

sub print_summary {
  my ($self, $prop) = @_;
  my $report = $prop->accession."\t".$prop->name."\t".$prop->result."\n";
  
  my $fh = $self->summaryFH;
  print $fh $report;
    
}

sub print_table {
  my ($self, $prop ) = @_; 
#acc     description     property_value  step_number     step_name       step_value      required?       evidence_type   evidence_name   best_HMM_hit    HMM_hit_count
#GenProp0001     chorismate biosynthesis via shikimate   YES     5       shikimate kinase        yes     required        HMM     PF01202 tr|F4HBD8|F4HBD8_GALAU  2
  my ($row, $report); 
  $row = join("\t", $prop->accession, $prop->name, $prop->result);
  foreach my $step (sort { $a->order cmp $b->order} @{ $prop->get_steps }){
    next if($step->skip);
    my $evRef = $step->get_evidence;
    my %et; #Store evidence type
    
    foreach my $e (@{ $step->get_evidence }){
      if($e->gp){
        $et{"Genome Property"}++;
      }else{
        $et{"InterPro"}++;
      }
    }
    my $evidence_type = join(",", keys(%et));
    
    $report .= $row."\t".$step->order."\t".$step->step_name."\t".$step->found."\t".(defined($step->required) ? "YES" : "NO")."\t$evidence_type\n";
  }
  #Now write the report to the filehandle.
  my $tfh = $self->tableFH;
  print $tfh $report;
  
}

sub print_compressed_table {

  my ($self, $prop ) = @_; 
  #acc     description     property_value [array of step values]
  #GenProp0001     chorismate biosynthesis via shikimate   [0,1,1,1,1] 
  my ($row, $report); 
  $row = join("\t", $prop->accession, $prop->name, $prop->result);
  my @stepRes;

  foreach my $step (sort { $a->order cmp $b->order} @{ $prop->get_steps }){
    if($step->skip){
      push(@stepRes, -1);    
    }else{
      push(@stepRes, $step->found);
    }
  }
  
  $row .= "\t[".join(",", @stepRes)."]";
  #Now write the report to the filehandle.
  my $tfh = $self->compTableFH;
  print $tfh $row;
  
}

sub print_long {
  my($self, $prop) = @_;

  my $report = "PROPERTY: ".$prop->accession."\n";
  $report .= $prop->name."\n";
  #TODO switch sort to <=>, once orders are replaced.
  foreach my $step (sort { $a->order cmp $b->order} @{ $prop->get_steps }){
    #	print LONGFORM (".\tSTEP NUMBER: $steps{$step_p}[2]\n.\tSTEP NAME: $steps{$step_p}[3]\n");
    $report .= ".\tSTEP NUMBER: ".$step->order."\n";
    $report .= ".\tSTEP NAME: ".$step->step_name."\n";  
    $report .= ".\t.\trequired\n" if (defined($step->required) and $step->required == 1);
    
    foreach my $ev (@{$step->get_evidence}){
      if($ev->interpro){
        $report .= ".\t.\tINTERPRO: ".$ev->interpro."; ".$ev->accession;
        $report .= " (clustered)" if($ev->clustered);
        $report .= "\n";
      }elsif($ev->gp){
        $report .= ".\t.\tGENPROP: ".$ev->gp."\n";
      }else{
        warn "Unknown step type\n";
        p($ev);
      }
    }
    #TODO: Does this relate to the step or evidence
    $report .= ".\tSTEP RESULT: ".($step->found == 1 ? 'yes' : 'no')."\n";
  
  }
    
  $report .= "RESULT: ".$prop->result."\n";
  my $fh = $self->longFH;  
  print $fh $report;

}	

sub print_matches {
  my($self, $prop) = @_;

  #Unlike the other views, which are genome property orientated,
  #this one is gene orientated.

  if($prop->result eq "YES" or $prop->result eq "PARTIAL"){
    my $fh = $self->matchFH;  
    my $pDESC = $prop->accession."\t".$prop->name;
    
    #TODO switch sort to <=>, once orders are replaced.
    foreach my $step (sort { $a->order cmp $b->order} @{ $prop->get_steps }){
      if($step->found == 1){
        my $sDESC = $step->order."\t".$step->step_name."\t";
        $sDESC .= "required" if (defined($step->required) and $step->required == 1);
          
        foreach my $ev (@{$step->get_evidence}){
          my $eDESC;
          my $seenSeq;
          if($ev->interpro){
            #Only report the clustered matches that passed the proximity check.
            my $seqsRef = $self->_uses_clustered_rule($ev) ? $step->clustered_matches($ev->accession)
                                                          : $self->get_family( $ev->accession );
            if($seqsRef){
              my @seqs = @{ $seqsRef };
              foreach my $s (@seqs){
                #next if($seenSeq->{$s});
                my $report = $ev->interpro."\t".$ev->accession."\t$s\t";
                print $fh $pDESC."\t$sDESC\t$report\n";
                #$seenSeq->{$s}++;
                #last;
                
              }
            }
          }
        }
      }
    } 
  }
}


sub open_outputfiles {
  my ($self) = @_;
    
   #LONGFORM_REPORT_1005058 TABLE_1005058 SUMMARY_FILE_1005058  
   if(!defined($self->{outdir})){
    $self->{outdir} = ".";  
   }
   
   my $root = $self->{outdir};
   
   #Open 
   foreach my $f (@{$self->{outfiles}}){
      warn "Opening filehandle based on $f\n";
      if($f eq 'long'){
        my $file = $root."/LONGFORM_REPORT_".$self->{name};
        my $fh;
        open($fh, '>', $file) or die "Failed to open $file:[$!]\n";
        $self->longFH($fh);
      }elsif($f eq 'summary'){
        my $file = $root."/SUMMARY_FILE_".$self->{name};
        my $fh;
        open($fh, '>', $file) or die "Failed to open $file:[$!]\n";
        $self->summaryFH($fh);
      }elsif($f eq 'table'){
        my $file = $root."/TABLE_".$self->{name};
        my $fh;
        open($fh, '>', $file) or die "Failed to open $file:[$!]\n";
        $self->tableFH($fh);
      }elsif($f eq 'comp_table'){
        my $file = $root."/COMP_TABLE_".$self->{name};
        my $fh;
        open($fh, '>', $file) or die "Failed to open $file:[$!]\n";
        $self->compTableFH($fh);
      }elsif($f eq 'web_json'){
        my $file = $root."/JSON_".$self->{name};
        my $fh;
        open($fh, '>', $file) or die "Failed to open $file:[$!]\n";
        $self->jsonFH($fh);
      }elsif($f eq 'match'){
        my $file = $root."/MATCHES_".$self->{name};
        my $fh;
        open($fh, '>', $file) or die "Failed to open $file:[$!]\n";
        $self->matchFH($fh);
      }else{
        die "Unknown output type $f\n";
      }
   }
   return 1;
}

sub longFH {
  my ( $self, $fh ) = @_;
  
  if($fh){
    if(ref($fh) eq "GLOB"){
      $self->{longFH} = $fh;    
    }else{
      croak("Filehandle not passed in\n"); 
    }
  }
  return ($self->{longFH});
}

sub summaryFH {
  my ( $self, $fh ) = @_;
  
  if($fh){
    if(ref($fh) eq "GLOB"){
      $self->{summaryFH} = $fh;    
    }else{
      croak("Filehandle not passed in\n"); 
    }
  }
  return ($self->{summaryFH});
}

sub tableFH {
  my ( $self, $fh ) = @_;
  
  if($fh){
    if(ref($fh) eq "GLOB"){
      $self->{tableFH} = $fh;    
    }else{
      croak("Filehandle not passed in\n"); 
    }
  }
  return ($self->{tableFH});
}

sub compTableFH {
  my ( $self, $fh ) = @_;
  
  if($fh){
    if(ref($fh) eq "GLOB"){
      $self->{compTableFH} = $fh;    
    }else{
      croak("Filehandle not passed in\n"); 
    }
  }
  return ($self->{compTableFH});
}

sub jsonFH {
  my ( $self, $fh ) = @_;
  
  if($fh){
    if(ref($fh) eq "GLOB"){
      $self->{jsonFH} = $fh;    
    }else{
      croak("Filehandle not passed in\n"); 
    }
  }
  return ($self->{jsonFH});
}

sub matchFH {
  my ( $self, $fh ) = @_;
  
  if($fh){
    if(ref($fh) eq "GLOB"){
      $self->{matchFH} = $fh;    
    }else{
      croak("Filehandle not passed in\n"); 
    }
  }
  return ($self->{matchFH});
}

sub removeSummaryFH {
  my ($self);
  delete $self->{summaryFH};
}

sub removeCompTableFH {
  my ($self);
  delete $self->{compTableFH};
}

sub removeTableFH {
  my ($self);
  $self->{tableFH}=undef;
}

sub removeLongFH {
  my ($self);
  $self->{longFH}=undef;
}

sub close_outputfiles {
  my($self) = @_;
  close($self->longFH) and $self->removeLongFH if ($self->longFH);
  close($self->summaryFH) if ($self->summaryFH); 
  $self->removeSummaryFH;
  close($self->tableFH) and $self->removeTableFH if ($self->tableFH);
  close($self->compTableFH) if($self->compTableFH);
  close($self->jsonFH) if($self->jsonFH);
  close($self->matchFH) if($self->matchFH);
  #TODO - remove these from the object.
  return(1);
}
	
sub debug{
  my($self, $bool) = @_; 
  $self->{debug} = $bool if($bool);
  return($self->{debug});
}

sub set_options {
  my($self, $options, $role) = @_;

  #Check the inputs are logical. 
  if(!$options->{seqs} and !$options->{matches} and 
      !$self->{seqs} and !$self->{matches} and $role eq 'cal'){
    croak("Need either a sequence file or match file"); 
  }

  my $input = 0;
  if (defined($options->{all})){
    $input++;
  }
  if ($options->{list}){
    $input++;
  }
  
  if ($options->{property}){
    $input++;
  }

  if ($input > 1){
	  die "Only one form of input please\n";
  }elsif($input == 0){
	  die "Please provide a form of input\n";
  }

  foreach my $k (keys %{$options}){
    $self->{$k} = $options->{$k};
  }

  if(defined($self->{clustered}) and $self->{clustered} ne 'strict' and $self->{clustered} ne 'permissive'){
    die "The clustered option must be strict or permissive\n";
  }
  if($self->{cluster_distance} !~ /^\d+$/){
    die "The cluster_distance option must be a whole number of base pairs\n";
  }

  #TODO: die on fatal option omitions.
  return 1;
}

sub _open_input_file {
  my ($self, $file) = @_;

  my $fh = IO::Uncompress::Gunzip->new($file, Transparent => 1, MultiStream => 1);
  croak("Failed to open $file: $GunzipError\n") unless $fh;

  return $fh;
}


sub define_sequence_set {
  my($self, $seq_blob) = @_;
  
  my @sb;
  if($seq_blob){
    @sb = split(/\n/, $seq_blob);
    foreach my $r (@sb){
      
    }
  }elsif($self->{seqs}){
    #open and read the sequence file
    my $fh = $self->_open_input_file($self->{seqs});
    @sb = <$fh>;
    close($fh);
  }


  foreach my $l (@sb){
    if($l =~ /^>(\S+)/){
      $self->{seqs_and_annotations}->{$1} = [];
    }
  }

}

sub get_sequence_set {
  my ($self) = @_;

  return($self->{seqs_and_annotations});
}

sub annotate_sequences {
  my($self) = @_;

  #Read from match file if defined
  if($self->signature_matches){
    ;
  }else{
    die "Can only run using pre-calculated i5 output\n";
  }
  $self->transform_annotations;
  $self->annotated(1);
}

sub signature_matches {
  my ($self) = @_;
  
  if(!$self->{read_sig} && $self->{matches}){
    my @mb;
    #If we get matches passed as a blob on command line via CGI,
    #split on carriage return, other than that, try and read
    #from the command line
    if(defined($self->{match_source}) and $self->{match_source} eq "inline"){
      @mb = split(/\n/, $self->{matches});  
    }else{
      my $fh = $self->_open_input_file($self->{matches});
      @mb = <$fh>;
      close($fh);
    }
    
    foreach my $line (@mb){
      chomp($line);
      my @i5 = split(/\t/, $line);
      push(@{ $self->{seqs_and_annotations}->{$i5[0]} }, $i5[4]);
    }
    $self->{read_sig} =1;  
  }
  return $self->{read_sig};
  
}

sub transform_annotations {
  my ($self) = @_; 

  my %fams;
  foreach my $seq (keys %{$self->{seqs_and_annotations} }){
    foreach my $fam (@{ $self->{seqs_and_annotations}->{$seq} }){
      push(@{$fams{$fam}}, $seq); 
    }
  }
  $self->families(\%fams);
}

sub families{
  my ($self, $fams) = @_;
  
  if($fams){
    $self->{families} = $fams;  
  }  
  
  return($self->{families});
}

sub get_family {
  my ($self, $acc) = @_;
  
  if(defined ($self->{families}) and defined ($self->{families}->{$acc})){
    return $self->{families}->{$acc};
    
  }else{
    return 0; 
  }
}

sub evaluate_properties {
  my ($self) = @_;  

  #GP can require other genome properties. We therefore need to
  #recursively search them, until all are evaluated.
  $self->incomplete(1);
  while ($self->incomplete) {
	  $self->check_results();
  }
  return 1;
}


sub incomplete {
  my($self, $value) = @_;
  
  if(defined($value)){
    $self->{incomplete} = $value;
  }
  
  return($self->{incomplete});
}

sub skip_def{
  my ($self, $type) = @_;
  
  if(!$type){
    croak("No property type passed in!\n");  
  }
  my $skip = 0;
  if(defined( $self->{_skip}->{$type} )){
    $skip = 1;  
  }
    
  return $skip;
}

sub check_results{
  my($self) = @_;
	
  my @prop_list = keys(% { $self->get_defs });  
  
  my $miss = 0;
  foreach my $acc (@prop_list){
    #Some property types are not required to be evaluated
    next if ( $self->skip_def( $self->get_def($acc)->type )); 
    
    #Evaluate the property, two scenarios, not tested, or tested and it had an unevaluated dependency
    if (!defined($self->get_def($acc)->evaluated) or $self->get_def($acc)->evaluated == 0 ){
		  $self->evaluate_property($acc);
		  
		  
		  #We should have evaluated the proprtery, unless it has an unevaluated deps.
      if ($self->get_def($acc)->evaluated == 0 ){
				$miss++;
			  if ($self->debug) {
          print "$acc is still missing\n";
			    print "value of missing: $miss\n";
        }
			}
		} 
	}
	
  if ($self->debug) {
    print "MISSING is $miss <-----------------------------------\n";
  }
  if($miss == 0){
    $self->incomplete(0);
  }
}


sub evaluate_property {
  my($self, $acc) = @_; 

	my $found = 0;          # count of REQUIRED steps that were found
	my $missing = 0;        # count of REQUIRED steps that were not found
	my $required_total = 0; # total number of REQUIRED steps (excluding skipped)
	my $all_found = 0;      # count of ALL steps (required + optional) that were found
	my $all_total = 0;      # count of ALL non-skipped steps

  if($self->debug){
    print "In evaluate proptery, $acc\n";
  }

  my $def = $self->get_def($acc);
  $def->evaluated(1);
  foreach my $step (@{ $def->get_steps }){
    if($self->debug){
      print "Working on step:".p($step)."\n";
    }

    $self->evaluate_step($step);

    #If we were not able to evaluate the step, bail
    if($step->evaluated == 0){
      $def->evaluated(0);
      return;
    }
  }

  #Clustered evidence depends on where the other steps were found,
  #so it can only be checked once every step has been evaluated.
  $self->evaluate_clustered_steps($def);

  foreach my $step (@{ $def->get_steps }){
    #Skipped steps (no evidence) are ignored entirely, matching the
    #pygenprop reimplementation which never creates a step without
    #functional elements/evidence.
    next if $step->skip == 1;

    $all_total++;
    if($step->found){
      $all_found++;
    }

    if($step->required){
      $required_total++;
      if($step->found){
        $found++;
      }else{
        $missing++;
      }
    }
  }


  #Assignment matches the pygenprop reimplementation:
  #  - If the property has one or more required steps, only the required
  #    steps are consulted (calculate_property_assignment_from_required_steps):
  #      YES     if every required step is found
  #      PARTIAL if the number of found required steps is greater than the threshold
  #      NO      otherwise
  #    Non-required (optional) steps do NOT contribute to the threshold
  #    comparison. This differs from the original Perl behaviour, which counted
  #    every found step toward the threshold and could therefore report PARTIAL
  #    based solely on optional-step evidence.
  #  - If the property has no required steps at all, every step is consulted
  #    (calculate_property_assignment_from_all_steps):
  #      YES     if every step is found
  #      NO      if no step is found
  #      PARTIAL otherwise
  if($required_total > 0){
    if($found == $required_total){
      $def->result('YES'); #All required steps found.
    }elsif($found > $def->threshold){
      $def->result('PARTIAL'); #Some required steps found, above threshold, but one or more missing.
    }else{
      $def->result('NO'); #At or below the threshold of found required steps.
    }
  }else{
    if($all_found == $all_total){
      $def->result('YES');
    }elsif($all_found == 0){
      $def->result('NO');
    }else{
      $def->result('PARTIAL');
    }
  }
}



sub evaluate_step {
  my($self, $step) = @_;
  my $some_evidence = 0; #Indicates whether the step is associated with some evidence, not really required.

  #We need to know if this step have been evaluated. If it has a dependency on the another GP, that
  #has not been tested at this point, we will set this to being 0. Note, for steps with no evidence,
  #this remains true, but there is nothing to evaluate against.
  $step->evaluated(1);
  if ($self->debug) {
      print "Evaluating ".$step->step_name."\n";
  }
  my $evRef = $step->get_evidence;

  #Each evidence is individually evaluated to a success/fail boolean. Evidence explicitly
  #flagged as "sufficient" (EV lines with a third ';'-delimited field, e.g.
  #"EV  IPR000001; PF00001; sufficient;") is considered strong enough evidence on its own to
  #prove a step. Evidence without that flag is not, on its own, considered strong enough.
  #
  #Assignment rule (matches the pygenprop reimplementation and resolves
  #https://github.com/ebi-pf-team/genome-properties/issues/30):
  #  - If the step has one or more sufficient evidences, the step succeeds if ANY sufficient
  #    evidence succeeds (non-sufficient evidences are not consulted in this case).
  #  - If the step has no sufficient evidences at all, the step succeeds only if ALL of its
  #    (non-sufficient) evidences succeed.
  my @sufficient_results;
  my @non_sufficient_results;

  if (!scalar(@{$evRef})){
    if($self->debug){
      print ("\tthis step has no evidence\n");
    }
    $step->skip(1);
  } else {

    EV: foreach my $evObj (@{$evRef}){
    #if (!$ev_id) {print RESULTS ("there is no ev_id\n"); next;} - Add to the object method
      if ($self->debug) {print ("\tevidence is ".p($evObj)."\n");}
      if ($self->debug) {print "\tev type is: ".$evObj->type."\n";}

      my $ev_succeed = 0;

      if($evObj->gp){
        if(defined($self->get_defs->{ $evObj->gp })){
          # For properties a PARTIAL or YES result is considered success
          if( $self->get_defs->{ $evObj->gp }->result eq 'YES' or
              $self->get_defs->{ $evObj->gp }->result eq 'PARTIAL' ){
              $ev_succeed = 1;
           }elsif($self->get_defs->{ $evObj->gp }->result eq 'UNTESTED'){
              $step->evaluated(0);
#Todo - need to check this bit. Some times a step can have two evidences, so need to check this is okay.
           }
        }
      }elsif($evObj->interpro){
        #Need to annotated the sequences
        if(!$self->annotated){
          $self->annotate_sequences
        }
        #Clustered evidence is handled by evaluate_clustered_steps, and takes
        #no part in the sufficient/non-sufficient rule below.
        next EV if($self->_uses_clustered_rule($evObj));
        #See if the accession has been found
        if($self->get_family( $evObj->accession ) ){
          $ev_succeed = 1;
        }
      }else{
        die "unknown evidence object type\n";
      }

      if(defined($evObj->type) and $evObj->type eq 'sufficient'){
        push(@sufficient_results, $ev_succeed);
      }else{
        push(@non_sufficient_results, $ev_succeed);
      }
    }
  }

  my $success = 0;
  if(@sufficient_results){
    #Any successful sufficient evidence is enough on its own.
    foreach my $result (@sufficient_results){
      if($result){
        $success = 1;
        last;
      }
    }
  }elsif(@non_sufficient_results){
    #With no sufficient evidence present, every (non-sufficient) evidence must succeed.
    $success = 1;
    foreach my $result (@non_sufficient_results){
      if(!$result){
        $success = 0;
        last;
      }
    }
  }

  if ($success){
    $step->found(1);
    #$step->matched_sequences($mSeqs);
    #Need to get the sequences in here.
  }
}

#Equivalent of the JCVI HMM-CLUST method: a gene matching clustered evidence
#only fills a step if it lies within cluster_distance of a gene that matched
#non-clustered evidence for a different step of the same property.
sub evaluate_clustered_steps {
  my($self, $def) = @_;

  return if($self->_clustered_mode eq 'permissive');

  my @clusteredSteps;
  foreach my $step (@{ $def->get_steps }){
    next if($step->found);
    if(grep { $_->interpro and $_->clustered } @{ $step->get_evidence }){
      push(@clusteredSteps, $step);
    }
  }
  return unless(@clusteredSteps);

  #Without gene positions clustered evidence cannot be tested, so it is ignored.
  my $positions = $self->gene_positions;
  return unless($positions);

  my %anchors;
  foreach my $step (@{ $def->get_steps }){
    foreach my $ev (@{ $step->get_evidence }){
      next if(!$ev->interpro or $ev->clustered);
      my $seqs = $self->get_family( $ev->accession ) or next;
      push(@{ $anchors{ $step->order } }, grep { $positions->{$_} } @{$seqs});
    }
  }

  foreach my $step (@clusteredSteps){
    my @others = map { @{ $anchors{$_} } } grep { $_ ne $step->order } keys %anchors;
    next unless(@others);
    foreach my $ev (grep { $_->interpro and $_->clustered } @{ $step->get_evidence }){
      my $seqs = $self->get_family( $ev->accession ) or next;
      SEQ: foreach my $s (@{$seqs}){
        next SEQ unless($positions->{$s});
        foreach my $anchor (@others){
          next if($anchor eq $s);
          if($self->_genes_within( $positions->{$s}, $positions->{$anchor} )){
            $step->add_clustered_match( $ev->accession, $s );
            $step->found(1);
            next SEQ;
          }
        }
      }
    }
  }
}

sub _uses_clustered_rule {
  my($self, $ev) = @_;
  return ($ev->clustered and $self->_clustered_mode eq 'strict') ? 1 : 0;
}

sub _clustered_mode {
  my($self) = @_;
  return $self->{clustered} if(defined($self->{clustered}));
  return $self->{gff} ? 'strict' : 'permissive';
}

sub _genes_within {
  my($self, $pos1, $pos2) = @_;

  return 0 if($pos1->{contig} ne $pos2->{contig});
  #Negative when the two genes overlap.
  my $gap = max($pos1->{start}, $pos2->{start}) - min($pos1->{end}, $pos2->{end});
  return ($gap <= $self->{cluster_distance}) ? 1 : 0;
}

#Reads CDS positions from a GFF3 file. Proteins are matched to the InterProScan
#sequence identifiers using the ID, Name, locus_tag or protein_id attributes.
sub gene_positions {
  my($self) = @_;

  return unless($self->{gff});
  if(!$self->{gene_positions}){
    my %positions;
    open(my $fh, '<', $self->{gff}) or croak("Could not open GFF file ".$self->{gff}."\n");
    while(my $l = <$fh>){
      last if($l =~ /^##FASTA/);
      next if($l =~ /^#/);
      chomp($l);
      my @col = split(/\t/, $l);
      next unless(scalar(@col) == 9 and $col[2] eq 'CDS');
      my $pos = { contig => $col[0], start => min($col[3], $col[4]), end => max($col[3], $col[4]) };
      foreach my $tag (qw(ID Name locus_tag protein_id)){
        if($col[8] =~ /(?:^|;)$tag=([^;]+)/){
          $positions{$1} = $pos;
        }
      }
    }
    close($fh);
    $self->{gene_positions} = \%positions;
  }
  return $self->{gene_positions};
}

sub read_properties {
  my($self) = @_;
  
  
  #We can read genome_properties from a variety of files
    #During the option initailisation we should have 
    #been told if we have a file or a database.
    #If not, we are got to fail!
  if($self->gp_source eq 'file'){
    $self->_read_properties_from_flatfile
  }elsif($self->gp_source eq 'db'){
    $self->_read_properties_from_db 
  }else{
    croak("No source for genome properties.");
  }
  #See if we have been told about an execution order for the GPs.
  #This is important as some GP rely in others.  Otherwise, we will
  #run iteratively. 
  $self->set_execution_order(); 
}



sub write_properties {
  my ($self) = @_;
  
  if($self->gp_out eq 'file'){
    $self->_write_properties_to_file
  }elsif($self->gp_out eq 'db'){
    $self->_write_properties_to_db 
  }else{
    croak("No source for genome properties.");
  }
  return 1;
}
=head1 set_execution_order

Title: 

=cut

sub set_execution_order {
  my ($self) = @_;

  if(!defined($self->{execution_order})){
    my $order;
    #If we have a db connection, use this
    #If we have an order file use that
    #Else, just use the sorted list of accessions.
    if($self->gp_source eq 'file' and $self->_order_file){
      $order = $self->_read_order_from_file;
    }elsif($self->gp_source eq 'db'){
      $order = $self->_read_order_from_db;
    }else{
      #No order, so we will have to rely on iterative execution.
      $order = [ sort{ $a cmp $b } keys(%{$self->{properties}}) ];
    }

    #set it in the object.
    $self->{execution_order} = $order;
  }

  return($self->{execution_order});
}


sub _order_file {
  my ($self) = @_;
  return($self->{eval_order});
}

sub _read_order_from_file {
  my ($self) = @_;

  my @order;
  open(FH, '<', $self->_order_file) or croak("Failed to open order file:".$self->_order_file."\n");
  while(<FH>){
    chomp;
    push(@order, $_);
  }
  close(FH);

  
  return(\@order);

}

sub _gp_flatfile {
  my ( $self ) = @_;
  return($self->{gpdir}.'/'.$self->{gpff});
}


sub _gp_def_file {
  my ($self) = @_;
  #now construct the full path.
  return($self->{gpdir}.'/'.$self->{_prop_def});
}


sub _gp_step_file {
  my ($self) = @_;
  #now construct the full path.
  return($self->{gpdir}.'/'.$self->{_prop_step});
}

sub _gp_step_ev_file {
  my ($self) = @_;
  #now construct the full path.
  return($self->{gpdir}.'/'.$self->{_step_ev_link});
}


sub _write_properties_to_file {
  my ($self) = @_;

  $self->_write_defs;
  $self->_write_steps;  
  $self->_write_step_evidences;
}

sub _read_properties_from_flatfile{
  my ( $self ) = @_;

  if( -e $self->_gp_flatfile){
   if( -s $self->_gp_flatfile ){ 
      GenomePropertiesIO::parseFlatfile($self, $self->_gp_flatfile);   
   }else{
      croak("The genome properties file has no size ".$self->_gp_flatfile."\n");
   }
  }else{
    croak("The geneome properties file, ".$self->_gp_flatfile." does not exist.\n");
  }
}


sub _read_properties_from_file {
  my ($self) = @_;
  
  #See if the file is present, if it is, open it and read it.
  my $fh;
  if(-e $self->_gp_def_file){
    if(-s $self->_gp_def_file){
      open($fh, $self->_gp_def_file) or croak("Failed to open ".$self->gp_def_file."\n");
      $self->_parse_defs($fh);
      close($fh)
    }else{
      croak("The genome properties file has no size ".$self->_gp_def_file."\n");
    }
  }else{
    croak("The geneome properties file, ".$self->_gp_def_file." does not exist.\n");
  }
  
  if(-e $self->_gp_step_file){
    if(-s $self->_gp_step_file){
      open($fh, $self->_gp_step_file) or croak("Failed to open ".$self->gp_step_file."\n");
      $self->_parse_step($fh);
      close($fh)
    }else{
      croak("The genome properties file has no size ".$self->_gp_step_file."\n");
    }
  }else{
    croak("The geneome properties file, ".$self->_gp_step_file." does not exist.\n");
  }


  if(-e $self->_gp_step_ev_file){
    if(-s $self->_gp_step_ev_file){
      open($fh, $self->_gp_step_ev_file) or croak("Failed to open ".$self->gp_step_ev_file."\n");
      $self->_parse_evidence($fh);
      close($fh)
    }else{
      croak("The genome properties file has no size ".$self->_gp_step_ev_file."\n");
    }
  }else{
    croak("The geneome properties file, ".$self->_gp_step_ev_file." does not exist.\n");
  }

}

sub _parse_defs {
  my ($self, $fh) = @_;


  while( my $line = <$fh>){
    chomp($line);
    #splits up the prop def file lines by tab 				
    # Line looks like this:
    # 1       chorismate biosynthesis via shikimate   PATHWAY GenProp0001     1       2       SIMPLE
	  my @split_line = split(/\t/, $line); 
    
    #Need some QC here to check that we at least have the correct data
    if(scalar(@split_line) != 7){
      $self->croak("Error parsing line $line.  Go ".scalar(@split_line).", yet expected 7.\n");
    }


    my $gpd = GenomeProperties::Definition->new( {'accession' => $split_line[3],
                           'type'      => $split_line[2],
                           'name'      => $split_line[1],
                           'map'       => $split_line[0],
                           'public'    => $split_line[4],
                           'threshold' => $split_line[5],
                           'method'    => $split_line[6]} );
    $self->add_def($gpd);
  }
}

sub fromDESC {
  my ($self, $data ) = @_;
  my $gpd = GenomeProperties::Definition->new( {'accession' => $data->{AC},
                           'type'      => $data->{TP},
                           'name'      => $data->{DE},
                           'map'       => '',
                           'public'    => 1,
                           'threshold' => $data->{TH},
                           'author'    => $data->{AU},
                           'method'    => 'SIMPLE',
                           'parents'   => $data->{PARENT},
                           'comment'   => $data->{CC}, #Previously call description.
                           'private'   => $data->{private}});
  if($data->{REFS}){
    $gpd->refs($data->{REFS});
  }
  if($data->{DBREFS}){
    $gpd->dbrefs($data->{DBREFS});
  }
  my $order = 0;
  my $seenSteps;
  
	foreach my $step (@{$data->{STEPS}}){
    $order++;
    if(defined( $seenSteps->{ $step->{SN} } )){
      die("Duplicate step number ".$step->{SN}."\n");
    }else{
      $seenSteps->{ $step->{SN} }++;
    }
    my $stepObj = GenomeProperties::Step->new({ order        => $step->{SN},
                                             step_name    => $step->{ID},
                                             required     => $step->{RQ},
                                             in_rule      => 1,
                                             alter_name   => $step->{DN} });
    $gpd->add_step($stepObj);
    foreach my $ev (@{$step->{EVID}}){
      my $evObj;
      if($ev->{gp}){
        $evObj = GenomeProperties::StepEvidence->new( {
                                                      genome_property => $ev->{gp}
                                                     } );
      }else{
        $evObj = GenomeProperties::StepEvidence->new( {
                                                      interpro   => $ev->{ipr},
                                                      accession  => $ev->{sig},
                                                      type       => $ev->{sc},
                                                      get_go     => $ev->{go}
                                                     } );

      }
      $stepObj->add_evidence($evObj);
    }
  }
  $self->add_def($gpd);
}

sub toDESC {
  my ($self, $path ) = @_;

  my $defs = $self->get_defs;
  
  my $descfile;
  if ($path) {
    $descfile = $path . "/DESC";
  }
  else {
    $descfile = "DESC";
  }
  open( D, ">$descfile" )
    or die "Could not open $descfile file for writing to\n";
  $Text::Wrap::unexpand = 0;
  $Text::Wrap::columns = 76;

  #binmode(D, ":utf8");

  while (my ($acc, $desc) = each %$defs) {
    print D wrap( "AC  ", "AC  ", $desc->accession );
    print D "\n";
    print D "DE  ".$desc->name."\n";
    print D wrap( "TP  ", "TP  ", $desc->type); 
    print D "\n";

    print D wrap( "AU  ", "AU  ", $desc->author);
    print D "\n";
    print D wrap( "TH  ", "TH  ", $desc->threshold);
    print D "\n";
    
    if(defined($desc->refs)){
      foreach my $ref ( @{ $desc->refs } ) {
        if ( $ref->{RC} ) {
          print D wrap( "RC   ", "RC   ", $ref->{RC} );
          print D "\n";
        }
        print D "RN  [" . $ref->{RN} . "]\n";
        print D "RM  " . $ref->{RM} . "\n";
        print D wrap( "RT  ", "RT  ", $ref->{RT} );
        print D "\n";
        print D wrap( "RA  ", "RA  ", $ref->{RA} );
        print D "\n";
        print D "RL  " . $ref->{RL} . "\n";
      }
    }
    if(defined($desc->dbrefs)){    
      foreach my $xref ( @{ $desc->dbrefs } ) {
        #Print out any comment
        if ( $xref->{db_comment} ) {
          print D wrap "DC  ", "DC  ", $xref->{db_comment};
          print D "\n";
        }
        if ( $xref->{other_params} ) {
          print D "DR  "
          . $xref->{db_id} . "; "
          . $xref->{db_link} . "; "
          . $xref->{other_params} . ";\n";
        } else {
          print D "DR  " . $xref->{db_id} . "; " . $xref->{db_link} . ";\n";
        }
      }
    }
    if(defined($desc->parents)){
      foreach my $parent (@{ $desc->parents}){
        print D "PN  $parent\n"; 
      }
    }
    print D wrap( "CC  ", "CC  ", $desc->comment);
    print D "\n" if($desc->comment);
    print D wrap( "**  ", "**  ", $desc->private);
    print D "\n" if($desc->private);

    foreach my $step (@{ $desc->get_steps }){
      print D "--\n";
      print D "SN  ".$step->order."\n";
      print D "ID  ".$step->step_name."\n";
      print D "DN  ".$step->alter_name."\n" if(defined($step->alter_name));
      print D "RQ  ".$step->required."\n";
      foreach my $ev (@{ $step->get_evidence}){
        print D "EV  ".$ev->interpro."; ".$ev->accession."; ".$ev->type.";\n";
        if($ev->get_go){
          foreach my $go (@{$ev->get_go}){
            print D  "TG  ".$go.";\n";
          }
        }
      }
    }
    print D "//\n"    
  }
}


sub toDESCHack {
  my ($self, $path ) = @_;

  my $defs = $self->get_defs;
 
  my %comments; 
  open(C, '<', '/tmp/gp.clean.tsv') or die;
  while(<C>){
    chomp;
    my($prop, $comment, $private) = split(/\t/, $_);
    $comments{$prop}->{c} = $comment;
    $comments{$prop}->{p} = $private;

  }

  open(I, '<', '/Users/rdf/Projects/InterPro/GenomeProperties/projects/GenomeProperties/data/flatfiles/my_dumps/gp2interpro.txt') or die;
  my %interpro;
  while(<I>){
    chomp;
    my ($a, $b, $sig, $ipr) = split(/\t/, $_);
    $interpro{$sig} = $ipr;
  }
  close(I);
  
  open(I, '<', '/Users/rdf/Projects/InterPro/GenomeProperties/projects/GenomeProperties/data/flatfiles/my_dumps/sufficient_components.tsv') or die;
  my %sufficient;
  while(<I>){
    chomp;
    my ($gp, $sig, $status) = split(/\t/, $_);
    $sufficient{"$gp$sig"} = $status if($status and $status =~ /\S+/);
  }
  close(I);

  open(I, '<', '/Users/rdf/Projects/InterPro/GenomeProperties/projects/GenomeProperties/data/flatfiles/my_dumps/getGo.tsv') or die;
  my $go;
  while(<I>){
    chomp;
    my ($gp, $goacc, $step) = split(/\t/, $_);
    $go->{$gp}->{$step} = $goacc;
  }
  close(I);
 
 #GenProp0113     IUBMB: Corrin Biosynthesis (part 3)     http://www.chem.qmul.ac.uk/iubmb/enzyme/reaction/tetrapyr/corrin3.html
  open(I, '<', '/Users/rdf/Projects/InterPro/GenomeProperties/projects/GenomeProperties/data/flatfiles/my_dumps/iubmb.tsv') or die;
  my $dbrefs;
  while(<I>){
    chomp;
    my ($gp, $dc, $p1, $p2) = $_ =~ /(GenProp\d+)\tIUBMB\: (.*?)\thttp.*?(\w+)\/(\w+)\.html/;
    push(@{$dbrefs->{$gp}}, { db_comment => $dc, db_id => 'IUBMB', db_link => $p1, other_params => $p2 });
  }
  close(I);
 #GenProp0001     Phenylalanine, Tyrosine and Tryptophan Biosynthesis     http://www.genome.ad.jp/kegg/pathway/map/map00400.html
  open(I, '<', '/Users/rdf/Projects/InterPro/GenomeProperties/projects/GenomeProperties/data/flatfiles/my_dumps/kegg2.tsv') or die;
  while(<I>){
    chomp;
    my ($gp, $dc, $p1) = $_ =~ /(GenProp\d+)\t(.*?)\thttp.*?(map\d+)/;
    if($gp){
      push(@{$dbrefs->{$gp}}, { db_comment => $dc, db_id => 'KEGG', db_link => $p1});
    }else{
      warn "$_\n";
    }
  }
  close(I);
 
  open(I, '<', '/Users/rdf/Projects/InterPro/GenomeProperties/projects/GenomeProperties/data/flatfiles/my_dumps/metacycall') or die;
  while(<I>){
    chomp;
    my ($gp, $dc, $p1) = $_ =~ /(GenProp\d+)\t(.*?)\thttp.*?object=(.*)/;
    if($gp){
      push(@{$dbrefs->{$gp}}, { db_comment => $dc, db_id => 'MetaCyc', db_link => $p1});
    }else{
      warn "$_\n";
    }
  } 
  
  open(I, '<', '/Users/rdf/Projects/InterPro/GenomeProperties/projects/GenomeProperties/data/flatfiles/my_dumps/Wikipedia') or die;
  while(<I>){
    chomp;
    my ($gp, $link ) = $_ =~ /(GenProp\d+)\tWikipedia\:(\S+)/;
    if($gp){
      push(@{$dbrefs->{$gp}}, { db_id => 'Wikipedia', db_link => $link});
    }else{
      warn "$_\n";
    }
  } 

  my $refs;
  my %refCount;
  open(I, '<', '/Users/rdf/Projects/InterPro/GenomeProperties/projects/GenomeProperties/data/flatfiles/my_dumps/gp.lit_refs') or die;
  while(<I>){
    chomp;
    my ($gp, $pmid ) = $_ =~ /(GenProp\d+)\tPMID\:(\d+)/;
    if($gp){

    my $url = "http://www.ebi.ac.uk/europepmc/webservices/rest/search?format=json&query=ext_id:$pmid";
    my $content = get($url);

    my $ref = from_json($content);
      foreach my $refObj (@{$ref->{resultList}->{result}}){
        next unless $refObj->{pmid} == $pmid;
        $refObj->{authorString} =~ s/(\.)$//;
        $refCount{$gp} = 0 if(!$refCount{$gp});
        $refCount{$gp}++;
        push(@{$refs->{$gp}}, { RN => $refCount{$gp}, 
                              RM => $pmid,
                              RT => $refObj->{title},
                              RA => $refObj->{authorString}.";",
                              RL => $refObj->{journalTitle}.". ".$refObj->{pubYear}.";".$refObj->{journalVolume}.":".$refObj->{pageInfo}."."
                              });
        last;
      }
    }else{
      warn "$_\n";
    }
  } 


  foreach my $acc (sort{ $a cmp $b } keys %$defs) {
    my $descfile;
    if ($path) {
      $descfile = $path . "/DESC";
    }
    else {
      $descfile = "DESC";
    }
    if(!-d "all/".$acc){
      mkdir("all/$acc");
    }
  open( D, ">all/$acc/$descfile" )
    or die "Could not open $descfile file for writing to\n";
  $Text::Wrap::unexpand = 0;
  $Text::Wrap::columns = 76;
 

    my $desc = $defs->{$acc};
    print D wrap( "AC  ", "AC  ", $desc->accession );
    print D "\n";
    print D wrap( "DE  ", "DE  ", $desc->name);
    print D "\n";
    print D wrap( "TP  ", "TP  ", $desc->type); 
    print D "\n";

    print D wrap( "AU  ", "AU  ", "Haft D");

    print D "\n";
    print D wrap( "TH  ", "TH  ", defined($desc->threshold));
    print D "\n";
    
    if($refs->{$desc->accession}){
      $desc->refs($refs->{$desc->accession});
    }
    if(defined($desc->refs)){
      foreach my $ref ( @{ $desc->refs } ) {
        if ( $ref->{RC} ) {
          print D wrap( "RC   ", "RC   ", $ref->{RC} );
          print D "\n";
        }
        print D "RN  [" . $ref->{RN} . "]\n";
        print D "RM  " . $ref->{RM} . "\n";
        print D wrap( "RT  ", "RT  ", $ref->{RT} );
        print D "\n";
        print D wrap( "RA  ", "RA  ", $ref->{RA} );
        print D "\n";
        print D "RL  " . $ref->{RL} . "\n";
      }
    }
    if($dbrefs->{$desc->accession}){
      $desc->dbrefs($dbrefs->{$desc->accession});
    }
    if(defined($desc->dbrefs)){    
      foreach my $xref ( @{ $desc->dbrefs } ) {
        #Print out any comment
        if ( $xref->{db_comment} ) {
          print D wrap "DC  ", "DC  ", $xref->{db_comment};
          print D "\n";
        }
        if ( $xref->{other_params} ) {
          print D "DR  "
          . $xref->{db_id} . "; "
          . $xref->{db_link} . "; "
          . $xref->{other_params} . ";\n";
        } else {
          print D "DR  " . $xref->{db_id} . "; " . $xref->{db_link} . ";\n";
        }
      }
    }
    if(defined($desc->parents)){
      foreach my $parent (@{ $desc->parents}){
        print D "PN  $parent\n"; 
      }
    }


    
    $desc->{comment} = $comments{$desc->accession}->{c};
    $desc->{private} = $comments{$desc->accession}->{p} if ($comments{$desc->accession}->{p} =~ /\S+/);
    
    if($desc->comment and $desc->comment =~ /\S+/){
      print D wrap( "CC  ", "CC  ", $desc->comment);
      print D "\n";
    }

    if($desc->private and $desc->private =~ /\S+/){
      print D wrap( "**  ", "**  ", $desc->private);
      print D "\n" 
    }
    
    foreach my $step (sort{ $a->order cmp $b->order} @{ $desc->get_steps }){
      print D "--\n";
      print D "SN  ".$step->order."\n";
      print D "ID  ".$step->step_name."\n";
      print D "DN  ".$step->alter_name."\n" if(defined($step->alter_name) and $step->alter_name ne  'NULL' and $step->alter_name ne 'DEFER');
      print D "RQ  ".$step->required."\n";
      foreach my $ev (@{ $step->get_evidence}){
        if($interpro{ $ev->accession }){
          $ev->interpro($interpro{ $ev->accession });
          $ev->type('wibble');
          print D "EV  ".$interpro{$ev->accession}."; ".$ev->accession.";";
          if($sufficient{$desc->accession.$ev->accession}){
            print D " sufficient;";
          }
          print D "\n";
        }elsif( $ev->accession =~ /GenProp/){
          print D "EV  ".$ev->accession.";\n";

        }else{
          #warn $ev->accession ." not found in mapping\n" if($ev->accession =~ /PF|TIGR/);
          print D "EV  IPRXXXX; ".$ev->accession.";";
          if($sufficient{$desc->accession.$ev->accession}){
            print D " sufficient;";
          }
          print D "\n";

        }
        #print D "EV  ".$ev->interpro."; ".$ev->accession."; ".$ev->type.";\n";
        if(defined($go->{$desc->accession}->{$ev->get_go})){
          print D  "TG  ".$go->{$desc->accession}->{$ev->get_go}.";\n";
        }
      }
    }
    print D "//\n"    
  }
}


sub _write_defs { 
  my ($self, $fh) = @_;
  
  foreach my $gp ( keys %{ $self->get_defs }){
    my $def = $self->get_def($gp);
    my $line = join("\t", $def->map, $def->name, $def->type, $def->accession, $def->public, $def->threshold, $def->method);
    print "$line\n";
  }
  
}

sub add_def {
  my ($self, $def) = @_; 
  # add the property to a hash, indexed by accession
  #add mapping key if present
  #map according to category
  $self->{properties}->{$def->{accession}} = $def;
  $self->{_defmap}->{$def->{map}}  = $def->{accession};
  $self->{_defcat}->{$def->{type}} = $def->{accession};
}

sub get_defs {
  my ($self) = @_;
  return($self->{properties});
}

sub get_def {
  my ($self, $acc) = @_;
  
  if(!$acc){
    carp("No accession passed in.\n"); 
  }
  
  if(!defined($self->{properties}->{$acc})){
    carp("No property found for this accession ($acc)\n");  
    return 0;
  }else{
    return($self->{properties}->{$acc});
  }
}


sub _parse_step {
  my($self, $fh) = @_;
  while (my $line = <$fh>){																			#go through step table line by line
    #76      1       5       shikimate kinase        1       1       shikimate kinase (EC  2.7.1.71)
    chomp($line); 
	  
    my @split_line = split(/\t/, $line);															#split up the line of step table by tab
    if(scalar(@split_line) != 7){
      croak("Unexpected number of elements in line $line, of step file\n");
    }

    #Now look at the file and build up.....
    
    if ($split_line[0] =~ /^\d/){		#is 2 number in each line of step table = step_link??????
			my $step = GenomeProperties::Step->new({ int_step_id  => $split_line[0],
                                               def_id       => $split_line[1],
                                               order        => $split_line[2],
                                               step_name    => $split_line[3],
                                               required     => $split_line[4],
                                               in_rule      => $split_line[5],
                                               alter_name   => $split_line[6]
                                               });
      $self->add_step($step);
		}
  }
}

sub _write_steps{
  my ($self) = @_;

  foreach my $gp ( keys %{ $self->get_defs }){
    my $def = $self->get_def($gp);
    next if(!defined($def));
    foreach my $step (sort {$a->order <=> $b->order}@{$def->get_steps}){
      my $line = join("\t", $step->int_step_id,
                            $step->def_id, 
                            $step->order, 
                            $step->step_name, 
                            $step->required, 
                            $step->in_rule,
                            $step->alter_name );
      print "$line\n";
    }
  }

}


sub add_step {
  my ( $self, $step ) = @_;
    
  #Look up the GP def that we need to attach this step to.
  my $acc = $self->{_defmap}->{$step->def_id};
  $step->def_acc($acc) if(!$step->def_acc);
  my $def = $self->{properties}->{$acc};

  $def->add_step($step);
  $self->{_stepmap}->{$step->int_step_id}=$acc;
}



sub _parse_evidence{
  my ($self, $fh ) = @_;



  while (my $line = <$fh>) {																			#go through it line by line
    chomp($line);
    my @split_line = split(/\t/, $line);															#split the line up by tabs
    my $ev = GenomeProperties::StepEvidence->new( { step_ev_id => $split_line[0],
                                                    step       => $split_line[1],
                                                    accession  => $split_line[2],
                                                    method       => $split_line[3],
                                                    get_go     => $split_line[4],
                                                      } );
      
    $self->add_step_evidence($ev);
  }
}

sub _write_step_evidences {
  my ( $self ) = @_;
  

  foreach my $gp ( keys %{ $self->get_defs }){
    my $def = $self->get_def($gp);
    next if(!defined($def));
    foreach my $step (sort {$a->order <=> $b->order}@{$def->get_steps}){
      foreach my $ev (@{$step->get_evidence}){
        my $line = join( "\t", $ev->step_ev_id, $ev->step, $ev->accession, $ev->type, $ev->get_go);  
        print "$line\n";
      }
    }  
  }


}


sub add_step_evidence {
  my ($self, $ev) = @_;
    
  my $s = undef;
  my $acc = $self->{_stepmap}->{$ev->{step}};
  my $def = $self->{properties}->{$acc};
  foreach my $step ( @{$def->get_steps} ){
    if($step->int_step_id == $ev->{step}){
      $step->add_evidence($ev);
      last;
    }
  }
  
}


sub gp_source {
  my ($self) = @_;

  my $s = '';
  if($self->{db}){
    $s = 'db'
  }else{
    $s = 'file';
  }

  return($s);
}

sub gp_out {
  my ($self) = @_;

  my $s = '';
  if($self->{db_out}){
    $s = 'db'
  }else{
    $s = 'file';
  }

  return($s);
}

sub checkConnectivity {
  my($self) = @_;

  my %c;
  my $top = $self->get_def("GenProp0065");
  $top->checkConnection($self, \%c);
  
  my $error = 0;
  foreach my $p (keys %{ $self->get_defs }){
    if(!$c{$p}){
      print STDERR "$p appears to be disconnected\n";
      $error++;
    } 
  }

  if($error){
    print STDERR "There are $error Genome Properties that are disconnected from the tree\n";
  }

}


1;
