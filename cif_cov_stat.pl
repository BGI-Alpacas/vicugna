#!/usr/bin/perl
use Getopt::Long;

sub usage{
    print <<USAGE;
Name:
    $0
Description:

Usage:
    perl $0 <list> <ref> <AA> <out> <out2> <out3>
Options:
	-L,--list	List of cif.gz file
	-R,--ref	A fasta file of spike protein
	-A,--AA		The abbreviation of Amino acid
	-S,--stat	The result of statistics
	-F,--file	The result of all chain
	-P,--pass	The result after filer
	help		print this help information
e.g
	perl $0 -L list -R P0DTC2.ref.fasta -A AA.ref.txt -S stat.xls -F all_file.xls -P pass_file.xls

USAGE
    exit 0;
}
my ($help,$list,$ref,$AA,$out,$out2,$out3);
GetOptions(
	"list|L:s" => \$list,
	"ref|R:s" => \$ref,
	"AA|A:s" => \$AA,
	"stat|S:s" => \$out,
	"file|F:s" => \$out2,
	"pass|P:s" => \$out3,
	"help|?"     => \$help
);
die &usage if (!defined $list || !defined $ref ||!defined $AA ||!defined $out ||!defined $out2 ||!defined $out3|| defined $help );

my $line;

###输入fa参考去换行符
open IN,"$ref";
while(<IN>)
	{
	chomp;
	next if ($_=~/^>/);
	$line.=$_;
	}
###氨基酸缩写对应
my %AB;
open IN1,"$AA";
while (<IN1>)
	{
	chomp;
	my @aa = split/\t/,$_;
	$AB{$aa[1]}=$aa[0];
	}
###参考文件位点与氨基酸对应
my (%has,%hash1,%hash2,%hash3);
my @tmp = split//,$line;
for (my $i=0;$i<=$#tmp;$i++)
	{
	$hash{$i+1}=$AB{$tmp[$i]};
	}

my (%count,%count1,%count2);
my @all_name;

###统计文件个数
my $b_filter;###过滤前文件个数
my $b_filter1;###过滤前的chain个数
my $a_f_3;###过滤有3个chain的长度大于800的文件个数
my $a_f_3_1;###过滤有3个chain的长度大于800的chain个数
my $a_f_300;###过滤掉变异超过300个的chain的文件个数
my $a_f_300_1;###过滤掉变异超过300个的chain个数

open IN2,"$list";
while(<IN2>)
	{
	chomp;
	my $name1 = (split/\//,$_)[-1];
	my $name = (split/\./,$name1)[0];###name为该cif名称
	$b_filter++;
	if ($_=~/\.gz$/)
		{
		open GZ,"gzip -dc $_|" or die $!;
		}
	else {
		open GZ,"$_" or die $!;
		}
	while (<GZ>)
		{
		chomp;
		my @bb = split/\s+/,$_;
		my $type;
		my $pos;
		if ($bb[0] eq "ATOM")###匹配ATOM的行
			{
			my $acid = $bb[5];###氨基酸
			if ($bb[4]=~/([\[a-zA-Z\]]+)([\d+]+)/)###格式整理，1000之后位点格式有问题，将位点分离出来
				{
				$type = $1;
				$pos = $2;
				}
			else {
				$type = $bb[6];
				$pos = $bb[8];
				}
			my $name_t = $name."_cha".$type;
			push (@all_name,$name_t);###所有样本名放入数组
			$pos=~s/[a-zA-Z]+//g;###位点格式再过滤去除可能携带的字母
			$hash1{$pos} = $pos;###确定样本存在的位点
			$hash2{$name_t}{$pos} = $acid;###确定样本中位点对应的氨基酸
			}
		}
	}
###样本名列表去重
my %ha;
my @uniq = grep {++$ha{$_} < 2} @all_name;

my @TP;
my (%com,%blank,%mis);###定义一致，空白，不一致
my (%LE,%LE1);###定义位点个数
my @all_name2;

open OUT,">$out";
foreach my $name (@uniq)
	{
	$com{$name} = 0;
	$mis{$name} = 0;
	$blank{$name} = 0;
	$b_filter1++;
	foreach my $pos (sort {$a <=> $b} keys %hash)
		{
		if (defined $hash2{$name}{$pos})
			{
			if($hash{$pos} ne $hash2{$name}{$pos})###统计样本中位点与参考位点氨基酸不一致的个数
				{
				$com{$name}++;
				$hash3{$name}.="\t$hash2{$name}{$pos}";###不一致位点输出样本氨基酸类型
				}
			else {
				$mis{$name}++;
				$hash3{$name}.="\t1";###位点一致输出结果为1
				}
			}
		else {
			$blank{$name}++;
			$hash3{$name}.="\t0";###位点没有参考中的位点输出0
			}
		}
	my $f_name = (split/_/,$name)[0];###切割出文件名
	$length = $com{$name} + $mis{$name};###文件长度
	if ($length > 800)
		{
		$LE1{$f_name}++;###统计每个文件中超过800的chain的个数
		push (@all_name2,$name);###过滤掉不超过800的chain
		}
	}

my @a_f_all;###过滤后的总表

my @a_F;
my @b_F;

foreach my $name (@all_name2)
	{
	my $f_name = (split/_/,$name)[0];###切割出文件名
	if ($LE1{$f_name} >= 3)###每个文件需有3个chain的长度大于800
		{
		$a_f_3_1++;
		push (@a_F,$f_name);###有3个chain的长度大于800的文件
		if ($com{$name} < 300)###过滤掉变异超过300个的chain
			{
			$a_f_300_1++;
			push (@b_F,$f_name);###过滤掉变异超过300个的文件
			push (@a_f_all,$name);###过滤后的总表
			my $type = (split/_/,$name)[1];
			push (@TP,$type);
			foreach my $pos (keys %hash1)
				{
				if (defined $hash2{$name}{$pos})
					{
					if($hash2{$name}{$pos} eq $hash{$pos})###统计一致的个数
						{
						$count{$pos}++;###统计各位点一致的总数
						$count2{$type}{$pos}++;###统计各位点各个chain一致的总数
						}
					}
				}
			}
		}
	}
#foreach my $cc (sort {$a <=> $b} keys %count)
open OUT2,">$out2";
open OUT3,">$out3";
for (my $i=0;$i<=$#tmp;$i++)
	{
	my $n = $i+1;
	print OUT "\t$n";###输出位点表头
	print OUT2 "\t$n";
	print OUT3 "\t$n";
	}
print OUT "\ntotal\t";
print OUT2 "\n";
print OUT3 "\n";

#foreach my $cc (sort {$a <=> $b} keys %count)
for (my $i=0;$i<=$#tmp;$i++)
	{
	my $n = $i+1;
	if(defined $count{$n})
		{
		print OUT "$count{$n}\t";###输出各位点统计的一致总数
		}
	else {
		print OUT "0\t";
		}
	}
print OUT "\n";

my %ha1;
my @uniq1 = grep {++$ha1{$_} < 2} @TP;

###输出各chain一致的总数
foreach my $tt  (@uniq1)
	{
	print OUT "$tt\t";
	#foreach my $cc (sort {$a <=> $b} keys %count)
	for (my $i=0;$i<=$#tmp;$i++)
		{
		my $n = $i+1;
		if(defined $count2{$tt}{$n})
			{
			print OUT "$count2{$tt}{$n}\t";
			}
		else {
			print OUT "0\t";
			}
		}
	print OUT "\n";
	}

###输出所有ID的结果
foreach (sort keys %hash3)
	{
	print OUT2 "$_$hash3{$_}\n";
	}
close OUT2;

foreach (sort @a_f_all)
	{
	print OUT3 "$_$hash3{$_}\n";
	}
close OUT3;

my %ha2;
my @uniq2 = grep {++$ha2{$_} < 2} @a_F;

my %ha3;
my @uniq3 = grep {++$ha3{$_} < 2} @b_F;

$a_f_3 = @uniq2;
$a_f_300 = @uniq3;

###屏幕输出统计的文件及chain数
print "total_file_num\tafter_filter_800\tafter_filter_300\n$b_filter\t$a_f_3\t$a_f_300\n$b_filter1\t$a_f_3_1\t$a_f_300_1\n";
