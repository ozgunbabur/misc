package org.panda.misc;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Scanner;
import java.util.Set;

/**
 * For submitting Mutex jobs to exacloud.
 * Created by babur on 2/8/2016.
 */
public class CondorSubmit
{
	final static boolean COMPUTE_FDR = true;
	static final String base = "/home/users/babur/mutex-data/CCLE";
	static final Set<String> avoid = new HashSet<>(Arrays.asList(""));

	static final String SUB_DIR = "/home/users/babur/cluster-files";
	public static void main(String[] args) throws IOException
	{
		submitRecursive(base);
	}

	public static void submitRecursive(String dir) throws IOException
	{
		if (COMPUTE_FDR) submitComplete(dir);
		else if (hasParametersFile(dir)) submitSimpleRun(dir);

		for (File f : new File(dir).listFiles())
		{
			if (f.isDirectory() && !avoid.contains(f.getName()))
				submitRecursive(f.getPath());
		}
	}

	private static boolean hasParametersFile(String dir)
	{
		File f = new File(dir);
		if (!f.exists() || !f.isDirectory()) return false;

		for (File file : f.listFiles())
		{
			if (file.getName().equals("parameters.txt")) return true;
		}
		return false;
	}

	public static void submitComplete(String dir) throws IOException
	{
		boolean hasInitialScores = false;
		boolean hasCompletedAlready = false;
		boolean hasParametersFile = false;
		int randCnt = 0;

		for (File file : new File(dir).listFiles())
		{
			if (file.getName().equals("ranked-groups.txt"))
			{
				hasInitialScores = true;
				Scanner sc = new Scanner(file);
				String line = sc.nextLine();
				sc.close();
				hasCompletedAlready = line.split("\t").length > 2;
			}
			else if (file.getName().equals("parameters.txt"))
			{
				hasParametersFile = true;
			}
			else if (file.getName().equals("randscores") && file.isDirectory())
			{
				for (File f : file.listFiles())
				{
					if (f.getName().endsWith(".txt")) randCnt++;
				}
			}
		}

		if (!hasParametersFile || hasCompletedAlready) return;

		System.out.print("Submitting " + dir + " ... ");
		String h = System.currentTimeMillis() + "";
		int i = 1;
		BufferedWriter writer;
		String A = null;
		String B = null;

		if (!hasInitialScores)
		{
			A = SUB_DIR + "/script-" + h + (i++) + ".sub";
			writer = new BufferedWriter(new FileWriter(A));
			writer.write(SCRIPT
				.replace("<DIR>", dir + " no-random-run")
				.replace("<CNT>", "1"));
			writer.close();
		}

		if (randCnt < 100)
		{
			B = SUB_DIR + "/script-" + h + (i++) + ".sub";
			writer = new BufferedWriter(new FileWriter(B));
			writer.write(SCRIPT
				.replace("<DIR>", dir + " random")
				.replace("<CNT>", "" + (100 - randCnt)));
			writer.close();
		}

		String C = SUB_DIR + "/script-" + h + (i++) + ".sub";
		writer = new BufferedWriter(new FileWriter(C));
		writer.write(SCRIPT
			.replace("<DIR>", dir)
			.replace("<CNT>", "1"));
		writer.close();

		if (A == null && B == null)
		{
			Runtime.getRuntime().exec("condor_submit " + C);
		}
		else
		{
			String D = SUB_DIR + "/script-" + h + i + ".sub";
			writer = new BufferedWriter(new FileWriter(D));
			if (A != null) writer.write("Job A " + A + "\n");
			if (B != null) writer.write("Job B " + B + "\n");
			writer.write("Job C " + C + "\n");

			if (A != null)
			{
				writer.write("Parent A ");
				if (B != null) writer.write(" Child B\n");
				else writer.write(" Child C\n");
			}
			if (B != null) writer.write("Parent B Child C");

			writer.close();

			Runtime.getRuntime().exec("condor_submit_dag " + D);
		}

		pause(100);
		System.out.println("ok");
	}

	public static void submitSimpleRun(String dir) throws IOException
	{
		System.out.print("Submitting " + dir + " ... ");
		String file = SUB_DIR + "/script-" + System.currentTimeMillis() + ".sub";
		BufferedWriter writer = new BufferedWriter(new FileWriter(file));
		writer.write(SCRIPT
			.replace("<DIR>", dir)
			.replace("<CNT>", "1"));
		writer.close();
		Runtime.getRuntime().exec("condor_submit " + file);
		pause(100);
		System.out.println("ok");
	}

	private static void deleteTempFiles()
	{
		for (File file : new File(SUB_DIR).listFiles())
		{
			if (file.getName().endsWith(".sub"))
			{
				file.deleteOnExit();
			}
		}
	}

	/**
	 * @param time in miliseconds
	 */
	public synchronized static void pause(long time)
	{
		Waiter w = new Waiter();
		w.bekle(time);
	}

	static class Waiter
	{
		private synchronized void bekle(long time)
		{
			try
			{
				this.wait(time);
			}
			catch (InterruptedException e)
			{
				e.printStackTrace();
			}
		}
	}

	public static final String SCRIPT = "dir=/home/users/babur/\n" +
		"\n" +
		"getenv = True\n" +
		"\n" +
		"#Program\n" +
		"executable=/usr/bin/java\n" +
		"\n" +
		"#Arguments to the program\n" +
		"arguments=-jar /home/users/babur/Projects/mutex/target/mutex.jar <DIR>\n" +
		"\n" +
		"#stdout\n" +
		"output=$(dir)stdout.sum\n" +
		"\n" +
		"#stderr\n" +
		"error=$(dir)stderr.sum\n" +
		"\n" +
		"#Condor log file\n" +
		"log=$(dir)log.log\n" +
		"\n" +
		"#Grab 5GB memory\n" +
		"request_memory=5120\n" +
		"\n" +
		"#Queue the job\n" +
		"queue <CNT>\n";
}
