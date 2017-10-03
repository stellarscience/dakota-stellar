package TestSuite;

import java.io.*;

/*
 * C/C++ Users Journal Sept 2000 <br>
 * The Simplest Automated Unit Test Framework That Could Possibly Work <br>
 * Chuck Allison <br>
 */
public abstract class Test
{
    int nPass = 0;
    int nFail = 0;
    BufferedWriter sink;

    public final void setSink(BufferedWriter sink)
    {
        this.sink = sink;
    }
    public final Writer getSink()
    {
        return sink;
    }
    public final void test(String label,
                           boolean condition)
        throws IOException
    {
        if (!condition)
            fail(label);
        else
            succeed();
    }

    public final void fail(String label)
        throws IOException
    {
        if (sink != null)
        {
            sink.write(getClass().getName() +
                       " failure: " + label);
            sink.newLine();
        }
        ++nFail;
    }

    public final void succeed()
    {
        ++nPass;
    }
    public final void flush()
        throws IOException
    {
        if (sink != null)
        {
            sink.flush();
        }
    }
    public final int report()
        throws IOException
    {
        if (sink != null)
        {
            sink.write("Test \"" +
                       getClass().getName() +
                       "\":");
            sink.newLine();
            sink.write("\tPassed: " + nPass);
            sink.newLine();
            sink.write("\tFailed: " + nFail);
            sink.newLine();
        }
        return nFail;
    }
    public final int getNumPassed()
    {
        return nPass;
    }
    public final int getNumFailed()
    {
        return nFail;
    }
    public final void close()
        throws IOException
    {
        sink.close();
    }
    protected void reset()
    {
        nFail = nPass = 0;
    }

    public abstract void run()
        throws IOException;
}

