#!/usr/bin/env perl6

my @function_file = (dir).map({ $_ if $_ ~~ /\.R/ or $_ ~~ /\.sh/ });
for @function_file -> $function_file {
    next if $function_file ~~ /^\./;
    my $m = $function_file.slurp ~~ /(\#.*?\n\s*\n)/;

    say $function_file.Str;
    say "```";
    $m[0].Str.lines.map( -> $l { say $l if $l ~~ /^\#/});
    say "```";
}
