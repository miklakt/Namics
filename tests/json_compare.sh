#!/usr/bin/env bash

compare_json_profiles() {
  local left="$1"
  local right="$2"
  local coord_tol="$3"
  local value_tol="$4"

  perl -MJSON::PP -e '
use strict;
use warnings;

sub read_json {
  my ($path) = @_;
  open my $fh, "<", $path or die "ERROR: cannot open $path: $!\n";
  local $/ = undef;
  my $text = <$fh>;
  close $fh;
  return JSON::PP::decode_json($text);
}

sub as_problem_object {
  my ($data, $path) = @_;
  if (ref($data) eq "ARRAY") {
    die "ERROR: JSON array is empty: $path\n" if scalar(@$data) == 0;
    $data = $data->[-1];
  }
  die "ERROR: JSON root is not an object: $path\n" unless ref($data) eq "HASH";
  if (exists $data->{problems}) {
    my $problems = $data->{problems};
    die "ERROR: JSON problems section is malformed: $path\n"
      unless ref($problems) eq "ARRAY" && scalar(@$problems) > 0;
    $data = $problems->[-1];
    die "ERROR: JSON problem entry is not an object: $path\n" unless ref($data) eq "HASH";
  }
  return $data;
}

sub extract_numeric_arrays {
  my ($obj) = @_;
  my %arrays = ();
  while (my ($key, $value) = each %$obj) {
    next unless ref($value) eq "ARRAY";
    next if scalar(@$value) == 0;
    my @numbers = ();
    my $ok = 1;
    for my $item (@$value) {
      if (ref($item)) { $ok = 0; last; }
      if ($item =~ /^-?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/) {
        push @numbers, 0.0 + $item;
      } else {
        $ok = 0;
        last;
      }
    }
    $arrays{$key} = \@numbers if $ok;
  }
  return \%arrays;
}

sub compare_arrays {
  my ($left_name, $right_name, $left, $right, $coord_tol, $value_tol) = @_;

  for my $key (keys %$left) {
    die "ERROR: output is missing column: $key\n" unless exists $right->{$key};
    my $lv = $left->{$key};
    my $rv = $right->{$key};
    my $left_n = scalar(@$lv);
    my $right_n = scalar(@$rv);
    die "ERROR: row count mismatch for column $key ($left_name: $left_n, $right_name: $right_n)\n"
      unless $left_n == $right_n;
    my $tol = ($key eq "x" || $key eq "y" || $key eq "z") ? $coord_tol : $value_tol;
    for (my $i = 0; $i < $left_n; ++$i) {
      my $delta = abs($lv->[$i] - $rv->[$i]);
      if ($delta > $tol) {
        my $row = $i + 1;
        die "ERROR: mismatch at row $row, column $key (tol=$tol) left=$lv->[$i] right=$rv->[$i]\n";
      }
    }
  }

  for my $key (keys %$right) {
    die "ERROR: reference is missing column: $key\n" unless exists $left->{$key};
  }
}

my ($left_path, $right_path, $coord_tol, $value_tol) = @ARGV;
die "ERROR: usage compare_json_profiles <left> <right> <coord_tol> <value_tol>\n"
  unless defined $value_tol;

my $left_data = as_problem_object(read_json($left_path), $left_path);
my $right_data = as_problem_object(read_json($right_path), $right_path);
my $left_arrays = extract_numeric_arrays($left_data);
my $right_arrays = extract_numeric_arrays($right_data);

compare_arrays($left_path, $right_path, $left_arrays, $right_arrays, 0.0 + $coord_tol, 0.0 + $value_tol);
' "$left" "$right" "$coord_tol" "$value_tol"
}
