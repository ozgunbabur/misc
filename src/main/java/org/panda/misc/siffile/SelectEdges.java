package org.panda.misc.siffile;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;

/**
 * For printing certain lines in a SIF file.
 */
public class SelectEdges
{
	public static void main(String[] args) throws IOException
	{
//		print("ITGB3", "SRC", "/home/ozgun/Analyses/CausalPath-paper/causal-priors.txt");
		printWithSource("PAG1", "/home/ozgun/Analyses/CausalPath-paper/causal-priors.txt");
	}

	public static void print(String source, String target, String file) throws IOException
	{
		Files.lines(Paths.get(file)).forEach(l ->
		{
			String[] t = l.split("\t");
			if (t[0].equals(source) && t[2].equals(target))
			{
				System.out.println(l);
			}
		});
	}
	public static void printWithSource(String source, String file) throws IOException
	{
		Files.lines(Paths.get(file)).forEach(l ->
		{
			String[] t = l.split("\t");
			if (t[0].equals(source))
			{
				System.out.println(l);
			}
		});
	}
}
