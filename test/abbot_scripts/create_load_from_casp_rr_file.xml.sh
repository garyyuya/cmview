#!/bin/sh


if [ -z "$1" ]
then
	echo "usage: $0 <casp rr file name>"
	exit 1
fi

caspfile=$1
basename=`basename $caspfile`
scriptfile=caspfiletest.$basename.xml

cat  <<ENDXML > $scriptfile
<?xml version="1.0" encoding="UTF-8"?>
<AWTTestScript desc="Load from Casp RR file">
  <component class="javax.swing.JMenuItem" id="CASP RR File..." index="6" parent="JPopupMenu Instance 2" text="CASP RR File..." window="Contact Map Viewer 6" />
  <component class="javax.swing.JButton" id="Cancel" index="3" parent="JPanel Instance 4" text="Cancel" window="Load from CASP RR file" />
  <component class="cmview.View" id="Contact Map Viewer" root="true" title="Contact Map Viewer" />
  <component class="cmview.View" id="Contact Map Viewer 2" root="true" title="Contact Map Viewer" />
  <component class="cmview.View" id="Contact Map Viewer 3" root="true" title="Contact Map Viewer" />
  <component class="cmview.View" id="Contact Map Viewer 4" root="true" title="Contact Map Viewer" />
  <component class="cmview.View" id="Contact Map Viewer 5" root="true" title="Contact Map Viewer" />
  <component class="cmview.View" id="Contact Map Viewer 6" root="true" title="Contact Map Viewer" />
  <component class="cmview.View" id="Contact Map of T0354" root="true" title="Contact Map of T0354" />
  <component class="javax.swing.JMenu" id="File" index="0" parent="JMenuBar Instance" text="File" window="Contact Map Viewer 4" />
  <component class="javax.swing.JLayeredPane" id="JLayeredPane Instance" index="1" parent="JRootPane Instance" window="Contact Map Viewer 2" />
  <component class="javax.swing.JLayeredPane" id="JLayeredPane Instance 2" index="1" parent="JRootPane Instance 2" window="Load from CASP RR file" />
  <component class="javax.swing.JMenuBar" id="JMenuBar Instance" index="3" parent="JLayeredPane Instance" window="Contact Map Viewer 3" />
  <component class="javax.swing.JPanel" id="JPanel Instance" index="0" parent="JLayeredPane Instance 2" window="Load from CASP RR file" />
  <component class="javax.swing.JPanel" id="JPanel Instance 2" index="0" parent="JPanel Instance" window="Load from CASP RR file" />
  <component class="javax.swing.JPanel" id="JPanel Instance 3" index="0" parent="JPanel Instance 2" window="Load from CASP RR file" />
  <component class="javax.swing.JPanel" id="JPanel Instance 4" index="1" parent="JPanel Instance" window="Load from CASP RR file" />
  <component class="javax.swing.JPopupMenu" id="JPopupMenu Instance" index="0" invoker="File" />
  <component class="javax.swing.JPopupMenu" id="JPopupMenu Instance 2" index="0" invoker="Load from" />
  <component class="javax.swing.JRootPane" id="JRootPane Instance" index="0" parent="Contact Map Viewer" />
  <component class="javax.swing.JRootPane" id="JRootPane Instance 2" index="0" parent="Load from CASP RR file" />
  <component class="javax.swing.JTextField" id="JTextField Instance" index="1" parent="JPanel Instance 3" window="Load from CASP RR file" />
  <component class="javax.swing.JMenu" id="Load from" index="1" parent="JPopupMenu Instance" text="Load from" window="Contact Map Viewer 5" />
  <component class="cmview.LoadDialog" id="Load from CASP RR file" parent="Contact Map Viewer 3" title="Load from CASP RR file" />
  <component class="javax.swing.JButton" id="Ok" index="1" parent="JPanel Instance 4" text="Ok" window="Load from CASP RR file" />
  <launch args="[]" class="cmview.Start" desc="CMView" method="main" />
  <sequence>
    <action args="CASP RR File..." method="actionSelectMenuItem" />
    <wait args="Load from CASP RR file" class="abbot.tester.ComponentTester" method="assertComponentShowing" />
    <action args="JTextField Instance,$caspfile" method="actionKeyString" />
    <action args="Ok" class="javax.swing.AbstractButton" method="actionClick" />
    <wait args="Contact Map of T0354" class="abbot.tester.ComponentTester" method="assertComponentShowing" />
  </sequence>
  <assert component="Contact Map of T0354" method="getTitle" value="/Contact Map of T0354.*/" />
  <terminate />
</AWTTestScript>
ENDXML
