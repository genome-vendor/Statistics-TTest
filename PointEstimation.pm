package Statistics::PointEstimation;
use strict;
use Carp;
use vars qw($VERSION @ISA $AUTOLOAD);
use Statistics::Distributions qw(chisqrdistr tdistr fdistr udistr uprob chisqrprob tprob fprob);
use Statistics::Descriptive;
use POSIX;


@ISA= qw (Statistics::Descriptive::Full);
$VERSION = '1.0';
my %confidence_interval=  #data related to confidence interval 
(

	"significance" => undef,
	"alpha" => undef,
	"df" =>undef,
	"standard_error" => undef,
	"t_value" =>undef, 
	"t_statistic" =>undef,
	"t_prob" =>undef,
	"delta" =>undef,
	"upper_clm" => undef,
	"lower_clm" =>undef,
	"valid"  =>undef
);


	
sub new{
	my $proto = shift;
	my $class = ref($proto) || $proto;
	my $self = $class->SUPER::new();  
	my %confidence=%confidence_interval;
	$self->{confidence}=\%confidence;
	bless ($self, $class);  
	return $self;
}

sub compute_confidence_interval{
	my $self=shift;
	croak "sample size must be >1 to compute the confidence interval \n" if($self->count()<=1);
	$self->{confidence}->{'significance'}=95 if (!defined($self->{confidence}->{'significance'}));
	$self->{confidence}->{df}=$self->count()-1;
	$self->{confidence}->{alpha}=(100-$self->{confidence}->{significance})/2;
	$self->{confidence}->{alpha}/=100;
	$self->{confidence}->{standard_error}=$self->standard_deviation()/sqrt($self->count());
	$self->{confidence}->{t_value}=abs tdistr($self->{confidence}->{df},$self->{confidence}->{alpha});
	$self->{confidence}->{delta}=$self->{confidence}->{t_value}*$self->{confidence}->{standard_error};

	$self->{confidence}->{upper_clm}=$self->mean() +$self->{confidence}->{delta};
	$self->{confidence}->{lower_clm}=$self->mean() -$self->{confidence}->{delta};
	$self->{confidence}->{t_statistic}=$self->{confidence}->{standard_error}
						?($self->mean()/$self->{confidence}->{standard_error}):0;
	$self->{confidence}->{t_prob}=1- abs (tprob($self->{confidence}->{df},-1*$self->{confidence}->{t_statistic})-tprob($self->{confidence}->{df},$self->{confidence}->{t_statistic})) ;
	$self->{confidence}->{valid}=1;
	return 1;

}
sub add_data{
	my $self = shift;
	my $aref;

	if (ref $_[0] eq 'ARRAY') {
		$aref = $_[0];
	}
	else {
		$aref = \@_;
	}
	my $significance=$self->{confidence}->{'significance'} if (defined($self->{confidence}->{'significance'}));
	$self->SUPER::add_data($aref);
	$self->{confidence}->{'significance'}=$significance;
	$self->compute_confidence_interval() if ($self->count()>1) ;

	return 1;

}
sub set_significance{   # set the significance level. usually 90, 95 or 99 
	my $self=shift;
	my $significance=shift;
	$self->{confidence}->{'significance'}=$significance if (($significance>0)&&($significance<100));
	$self->compute_confidence_interval() if($self->count()>1);
	return 1;

}

sub print_confidence_interval{
	my $self=shift;
	print "mean:",$self->mean(),"\n";
	print "variance:",$self->variance(),"\n";
	my $confidence=$self->{confidence};

	foreach my $k ( keys %$confidence)
	{
		print "$k: $confidence->{$k} \n";
	}
	return 1;

}

sub output_confidence_interval{
	my $self=shift;
	croak "sample size must be >1 to compute the confidence interval\n" if($self->{confidence}->{valid}!=1);
	my $title=shift;
	print "Summary  from the observed values of the sample $title:\n";
	print "\tsample size= ", $self->count()," , degree of freedom=", $self->df(), "\n";
	print "\tmean=", $self->mean()," , variance=", $self->variance(),"\n";
	print "\tstandard deviation=", $self->standard_deviation()," , standard error=", $self->standard_error(),"\n";
	print "\t the estimate of the mean is ", $self->mean()," +/- ",$self->delta(),"\n\t",
		" or (",$self->lower_clm()," to ",$self->upper_clm," ) with ",$self->significance," % of confidence\n"; 
	print "\t t-statistic=T=",$self->t_statistic()," , Prob >|T|=",$self->t_prob(),"\n";
}

sub AUTOLOAD{
	my $self = shift;
	my $type = ref($self)
	or croak "$self is not an object";
	my $name = $AUTOLOAD;
	$name =~ s/.*://;     
	return if $name eq "DESTROY";
	if (exists $self->{_permitted}->{$name} ) {
		return $self->{$name};
	}
	elsif(exists $self->{confidence}->{$name})
	{
		return $self->{confidence}->{$name};
	}
	else
	{
		croak "Can't access `$name' field in class $type";
	}
}
1;


__END__

=head1 NAME

Statistics::PointEstimation - Perl module for computing the confidence interval in parameter estimation with Student's T distribution

=head1 SYNOPSIS

  use Statistics::PointEstimation;

  my @r=();
  for($i=1;$i<=32;$i++) #generate a uniformly distributed sample with mean=5   
  {

	  $rand=rand(10);
	  push @r,$rand;
  }

  my $stat = new Statistics::PointEstimation;
  $stat->set_significance(95); #set the significance(confidence) level to 95%
  $stat->add_data(@r);
  $stat->output_confidence_interval(); #output summary
  $stat->print_confidence_interval();  #output the data hash related to confidence interval estimation

  #the following is the same as $stat->output_confidence_interval();
  print "Summary  from the observed values of the sample:\n";
  print "\tsample size= ", $stat->count()," , degree of freedom=", $stat->df(), "\n";
  print "\tmean=", $stat->mean()," , variance=", $stat->variance(),"\n";
  print "\tstandard deviation=", $stat->standard_deviation()," , standard error=", $stat->standard_error(),"\n";
  print "\t the estimate of the mean is ", $stat->mean()," +/- ",$stat->delta(),"\n\t",
  " or (",$stat->lower_clm()," to ",$stat->upper_clm," ) with ",$stat->significance," % of confidence\n";
  print "\t t-statistic=T=",$stat->t_statistic()," , Prob >|T|=",$stat->t_prob(),"\n";




=head1 DESCRIPTION

  This module is a subclass of Statistics::Descriptive::Full. It uses T-distribution for point estimation 
  assuming the data is normally distributed or the sample size is sufficiently large. It overrides the 
  add_data() method in Statistics::Descriptive to compute the confidence interval with the specified significance
   level (default is 95%). It also computes the t-statistic=T and Prob>|T| in case of hypothesis 
  testing of paired T-tests.
 

=head1 AUTHOR

Yun-Fang Juan , Yahoo! Inc.  (yunfang@yahoo-inc.com)

=head1 SEE ALSO

Statistics::Descriptive Statistics::Distributions

=cut 
