package TestSuite;

import java.io.*;
import java.util.*;

/*
 * C/C++ Users Journal Sept 2000 <br>
 * The Simplest Automated Unit Test Framework That Could Possibly Work <br>
 * Chuck Allison <br>
 */
public final class Suite
{
    String name;
    BufferedWriter sink;
    ArrayList tests = new ArrayList();

    public Suite(String name)
        throws IOException
    {
        this.sink =
            new BufferedWriter(
                new OutputStreamWriter(System.out)
                              );
        this.name = name;
    }
    public Suite(String name,
                 BufferedWriter sink)
    {
        this.sink = sink;
        this.name = name;
    }

    public void addTest(Test t)
    {
        t.setSink(sink);
        tests.add(t);
    }

    public void run()
        throws IOException
    {
        for (int i = 0; i < tests.size(); ++i)
        {
            Test t = (Test) tests.get(i);
            t.run();
        }
    }
    public void report()
        throws IOException
    {
        for (int i = 0; i < tests.size(); ++i)
        {
            Test t = (Test) tests.get(i);
            t.report();
        }
        sink.flush();
    }
    public void flush()
        throws IOException
    {
        sink.flush();
    }

    public void close()
        throws IOException
    {
        sink.close();
    }
}
