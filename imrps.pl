#!/usr/bin/perl

my $pdb_input = $ARGV[0];
chomp($pdb_input);
$pdb_input = lc($pdb_input);
my $pdb_file = "/data/pdb/pdb$pdb_input.ent";
#my $pdb_file = "/data/pdb/pdb$pdb_input.ent";
###accessing the pdb file 
open my $myfile,"<","$pdb_file";
open my $outfile,">","atom_fasta.txt";
open my $out,">","seqres_fasta.txt";
$mod=0;
while (<$myfile>)
{
	if(/^NUMMDL/)	######## multiple models
	{
		$mod++;
	}
	if (/^HEADER\s+(.+)\s*\d\d\-\w+\-\d+\s+(\w+)/)
	{
		$id = $2;
		print $id,"\n";	######### PDB ID
		$name = $1;$name =~ s/\'/\\'/g;$name =~ s/\s+/ /g;
		print $name,"\n";	######### name of the protein
	}
	if(/^SOURCE\s+.+ORGANISM_SCIENTIFIC: (.+)(;)\s+$/)
	{	
		$organism = $1;$organism =~ s/\'/\\'/;
		$organism =~ s/\s+/ /g;
		print $organism,"\n";
			########## organism name
	}
	if(/^SOURCE\s+.+ORGANISM_SCIENTIFIC: (.+)(;{0})\s+$/)
	{
		$organism = $1;$organism =~ s/\'/\\'/;
		$organism =~ s/;//;$organism =~ s/\s+/ /g;
			########## organism name
	}
	if(/^SOURCE\s+3\s+(.+)/)
	{
		$hgj = $1;
		if($hgj !~ /^\w+.+:/)
		{
			$organism .= "$hgj"; $organism =~ s/;//;$organism =~ s/\s+/ /g;
			print "$organism\n";
		}
	}
	if(/^EXPDTA\s+(.+)/)
	{
		$exp_method = $1; $exp_method=~ s/\s+/ /g;	########## experimental method
		print "$exp_method\n";
	}
	if(/^REMARK\s+3\s+RESOLUTION RANGE HIGH \(ANGSTROMS\) : (.+)/)#REMARK   3   RESOLUTION RANGE HIGH (ANGSTROMS) : 2.20
	{
		$resolution = $1; $resolution =~ s/\s+/ /g;	########### resolution
		print "$resolution\n\n";
	}
	if(/^REMARK\s+3\s+R VALUE\s+\(WORKING SET\)\s+:\s+(.+)/)
	{
		$r_factor = $1;	$r_factor=~ s/\s+/ /g;	########### resolution
		print "$r_factor\n\n";
	}
	if($r_factor !~ /\w+/)
	{
		if(/^REMARK\s+3\s+R VALUE\s+\(WORKING SET, NO CUTOFF\)\s+:\s+(.+)/)
		{
			$r_factor = $1;	$r_factor=~ s/\s+/ /g;	########### resolution
			print "$r_factor\n\n";
		} 
	}
	if(/^SEQRES\s+\d+\s+(\w+)\s+\d+\s+(.+)/)
 	{
		$res = $2;
		$ch1 = $1;
		@ab = split(/ /,$res);
		foreach $x(@ab)
		{
			$k=$x.'###';	
			if($ch1 ne $ch2 && $ch1 ne '')   #####comparing the chains
			{
				$pam .= "%%%>$id:$ch1|SEQRES!!!\n";
				$pam .= $k;
				$ch2 = $ch1;
			}
			else
			{
				$pam .= $k;   #####using '###' to separate each residue  
			}
		}
	}
	if(/^MODRES\s+\w+\s+(\w+\s+\w+\s+\d+\s+\w+)/)
	{
		#print $1,"\n";
		push(@mod_data01,$1);
	}
	if(/^HETNAM\s+/)
	{
		$het = substr($_,11);
		push(@het_arr, $het);
		#HETNAM     TYS O-SULFO-L-TYROSINE
		
	}
}
@het_atm = @het_arr;
$a=0;
for $ha($a..scalar(@het_arr)-1)
{
	$rf = substr($het_arr[$ha],0,3);
	$lf = substr($het_arr[$ha],4);
	$b=$a+1;
	for $ah($b..scalar(@het_arr)-1)
	{
		$fr = substr($het_arr[$ah],0,3);
		$fl = substr($het_arr[$ah],4);
		if($rf eq $fr)
		{
			chomp($lf);$lf=~s/\s+$//;
			chomp($fl);$fl=~s/\s+$//;
			$lf .= $fl;
			$fnn = "$rf $lf\n";
			splice(@het_atm,$a,2,$fnn);
			#print "$a $b\n";
		}
	}
	$a++;
}

@mod_data01 = sort{substr($a,4,1) cmp substr($b,4,1)}@mod_data01;

foreach $f(0..scalar(@mod_data01)-1)
{
	$mod_data01[$f]=~/\w+\s+\w+\s+\d+\s+(\w+)/;
	if(length($1)==3)
	{
		push(@mod_data,$mod_data01[$f]);
	}
	
}
foreach $jj(@mod_data)
{
	if($jj=~/\w+\s+(\w)\s+\d+\s+\w+/)
	{
		push(@mod_cha,$1);
	}
	$m_re = substr($jj,0,3);
	foreach $h(@het_atm)
	{
		$re_m = substr($h,0,3);
		$nam_mod = substr($h,4);chomp($nam_mod);$nam_mod=~s/\s+$//;
		if($m_re eq $re_m)
		{
			$jj .= " $nam_mod";
		}
	}
}

push(@mod_chain,grep{! $mc{$_}++} @mod_cha);

@cann = split(/%%%/,$pam);  	#####splitting chainwise
shift(@cann);
foreach $e(@cann)
{
	$tn = 0;
	@sn1 = split(/!!!/,$e);    	 #####splitting and assigning fasta header as one element and residues as other element of an array @sn1
	@sn = split(/###/,$sn1[1]);	#########splitting and assigning each resiude as single element of an array @sn
	$tot_res = @sn;
	foreach $y(@sn)
	{
		$jjj = length($y);
		if($jjj <= 2)   ######### getting rid of nucleic acid#########
		{
			$tn++;
		}
	}
	if($tn == 0)	#########amino acid only
	{	
		$s_seq .= "\n$sn1[0]\n";
		if($sn1[0]=~/^>\w+:(\w+)\|/) 		####### getting only the amino acid chains
		{
			$chnan = $1;
			$tot_res_chn = "$chnan $tot_res";
			push(@chainss,$tot_res_chn);
		}
	}
}
if($mod == 0) 	############## single model is present
{
	open my $infile1, "<", "$pdb_file";
	while (<$infile1>)
	{
		if(/^TER/)
		{
			$last = $_;   ###assigning the last TER to $last 
		}
		
	}
	close($infile1);

	open my $mfile,"<","$pdb_file";
	while (<$mfile>)
	{

		######### extracting the atom record ##########

		if (/^SCALE3.+/..(/$last/))   ###using the flip-flop operator accessing the contents between SCALE3 and last TER
		{
			if(/^ATOM|HETATM\s{0,5}(.+)/)
			{
				push(@unique,$_);###pushing the atom and inbetween hetatm record to array: @unique
			}
		}
	}
	close($mfile);
}
else 	########## multiple models are present
{
	open my $infile2, "<", "$pdb_file";
	while (<$infile2>)
	{
		if (/^MODEL\s+1\s+/..(/^MODEL\s+2\s+/))###using the flip-flop operator accessing the contents between SCALE3 and last TER
		{
			if(/^ATOM|HETATM\s{0,5}(.+)/)
			{
				push(@unique,$_);
			}
		}
	}
	close($infile2);
}

foreach $q(@unique)
{	
	$mim = substr($q,17,10);
	if($mim =~ /(\w+)\s{0,1}(\w)\s{0,4}(-\d+\w?|\d+\w?)/)
	{
		$ress = $1;
		$chain = $2;
		$res_no = $3;
		$seq = "$ress $chain $res_no";
		push(@uni1,$seq);
	}
}

########eliminating the duplicates########

push(@unique1,grep{! $nn{$_}++} @uni1);
$val1 = $val2 = "";
@unique1_1 = @unique1;
for $y(0..scalar(@unique1_1)-2)
{
	$unique1_1[$y] =~ /\w+\s+(\w)\s+(-\d+\w?|\d+\w?)/;
	$no_r1 = $1;
	$reno1 = $2;
	$unique1_1[$y+1] =~ /\w+\s+(\w)\s+(-\d+\w?|\d+\w?)/;
	$no_r2 = $1;
	$reno2 = $2;
	if($reno1 eq $reno2 && $no_r1 eq $no_r2)
	{
		splice(@unique1,$y,1);
	}
}

  #####concatenating the atom record and separate atom record based on chain using &&& 
foreach $r(@unique1)
{
	$r =~ /(\w+)\s+(\w)\s+(-\d+\w?|\d+\w?)/;
	$val1 = $2;
	$res = $1;
	$no = $3;
	if($no!~/\d+$/)
	{
		push(@insert,$r);
	}
	if($val1 ne $val2 && $val1 ne '')    ###comparing the chains
	{
		$num .= "%%%>$id:$val1|ATOM!!!\n";	######extracting the residue number
		$num .= $no .'###';
		$val2 = $val1;
	}
	else
	{
		$num .= $no .'###';
	}
}
foreach $in(@insert)
{
	if($in=~/(\w+)\s+(\w)\s+\w+/)
	{
		if(length($1) == 3)
		{	
			push(@ins_cha,$2);
		}
	}
}
@try = split(/%%%/,$num);
shift(@try);
foreach $re(@chainss)
{
	$ne= substr($re,0,1);
	$de= substr($re,2); #total number of resiudes in a chain
	foreach $r(@try)
	{
		if($r=~/^>\w+\:(\w+)\|.+\n(\w+)\#.+\#(\d+)\#{3}$/)
		{
			$chain_loc = "$1 $2 $3 $de";	#$1= chain; $2-$3= start and end number; $de= total number of resiudes 
			if($ne eq $1)
			{
				push(@pdb_detail, $chain_loc);
			}
			
		}
	}
}
push(@ins_chain,grep{! $ic{$_}++} @ins_cha);
push(@req_cha, @mod_chain);push(@req_cha, @ins_chain);
push(@req_chain,grep{! $rc{$_}++} @req_cha);

if(scalar(@mod_chain)!=0)
{
	print "PDB_ID: $id\n";
	print "MODIFIED\n";
	foreach $rq(@mod_chain)
	{
		foreach $l(@pdb_detail)
		{
			$l=~/(\w+)\s+(\w+)\s+(\w+)\s+(\w+)/;
			$loc = "$2-$3";
			if ($rq eq $1)
			{
				$loc = "$2-$3";
				$seq_len = $4;
				print "CHAIN: $1\nRESIDUE_NO: $loc\nTOTAL_RESIDUE: $4\n";
			}
		}
		@mod_residue=@mod_position=@mod_name='';
		$modf=0;
		foreach $m(@mod_data)
		{
			if($m=~/(\w+)\s+(\w)\s+(-\d+\w?|\d+\w?)\s+(\w+)\s+(.+)/)
			{
				if($rq eq $2)
				{
					$vall1 = $2;
					$ress = $1;
					$noo = $3;
					$ress2 = $4;
					$m_nam= $5;
					push(@mod_residue,"$ress=>$ress2");
					push(@mod_name,"$ress=>$m_nam;");
					push(@mod_position, $noo);
					$modf++;
				}
			}
			if($m=~/(\w+)\s+(\w)\s+(-\d+\w?|\d+\w?)\s+(\w+)$/)
			{
				if($rq eq $2)
				{
					$vall1 = $2;
					$ress = $1;
					$noo = $3;
					$ress2 = $4;
					$m_nam= $5;
					push(@mod_residue,"$ress=>$ress2");
					push(@mod_name,"-");
					push(@mod_position, $noo);
					$modf++;
				}
			}
		}
		shift(@mod_residue);shift(@mod_position);shift(@mod_name);
		print "RESIDUES: @mod_residue\n";
		print "POSITION: @mod_position\n";
		print "TOTAL_MODIFIED: $modf\n";
		print "MODIFIED_RESIDUE: @mod_name\n\n";
		
	}
}

print "\n#################################################################################################\n\n";

if(scalar(@ins_cha)!=0)
{
	print "PDB_ID: $id\n";
	print "INSERTED\n";
	foreach $rq(@ins_chain)
	{
		foreach $l(@pdb_detail)
		{
			$l=~/(\w+)\s+(\w+)\s+(\w+)\s+(\w+)/;
			$loc = "$2-$3";
			if ($rq eq $1)
			{
				print "CHAIN: $1\nRESIDUE_NO: $loc\nTOTAL_RESIDUE: $4\n";
			}
		}
		@ins_residue=@ins_position='';
		$ins=0;
		foreach $ir(@insert)
		{
			if($ir=~/(\w+)\s+(\w)\s+(-\d+\w|\d+\w)/)
			{
				if($rq eq $2)
				{
					$val1 = $2;
					$res = $1;
					$no = $3;
					push(@ins_residue,$res);
					push(@ins_position,$no);
					$ins++;
				}
			}
		}
		shift(@ins_residue);shift(@ins_position);
		print "RESIDUES: @ins_residue\n";
		print "POSITION: @ins_position\n";
		print "TOTAL_INSERTED: $ins\n\n";
	}
}


