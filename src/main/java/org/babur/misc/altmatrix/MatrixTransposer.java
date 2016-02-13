package org.babur.misc.altmatrix;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.Scanner;

/**
 * For transposing a tab-delimited text file.
 *
 * Created by babur on 2/8/2016.
 */
public class MatrixTransposer
{
	public static void run(String inFile, String outFile) throws IOException
	{
		Scanner sc = new Scanner(new File(inFile));

		List<String>[] row = null;

		while (sc.hasNextLine())
		{
			String[] token = sc.nextLine().split("\t");
			if (row == null)
			{
				row = new List[token.length];
				for (int i = 0; i < row.length; i++)
				{
					row[i] = new ArrayList<>();
				}
			}

			for (int i = 0; i < token.length; i++)
			{
				row[i].add(token[i]);
			}
		}
		sc.close();

		if (row == null)
		{
			System.err.println("No data in the input file");
			return;
		}

		BufferedWriter writer = new BufferedWriter(new FileWriter(outFile));

		for (int i = 0; i < row.length; i++)
		{
			if (i > 0) writer.write("\n");

			for (int j = 0; j < row[i].size(); j++)
			{
				writer.write(row[i].get(j));
				if (j < row[i].size() - 1) writer.write("\t");
			}
		}

		writer.close();
	}

	public static void main(String[] args) throws IOException
	{
		String dir = "C:/Users/babur/Documents/mutex/CCLE/";
//		String dir = "C:\\Users\\babur\\Documents\\mutex\\CCLE\\";
		run(dir + "CCLE_mutation_featurematrix.tab", dir + "DataMatrix.txt");
	}
}
