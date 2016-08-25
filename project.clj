(defproject clatrix/clatrix "0.5.1-SNAPSHOT"
  :description "Because using matrices in Clojure needs to not suck."
  :url "http://tel.github.com/clatrix"
  :license {:name "MIT License"
            :url "http://www.opensource.org/licenses/MIT"}
  :resource-paths ["native"]
  :plugins [[lein-expectations "0.0.8"]]

  :dependencies [[org.clojure/clojure "1.8.0"]
                 [slingshot "0.12.2"]
                 [org.jblas/jblas "1.2.3"]
                 [net.mikera/core.matrix "0.36.1"]]
  
  :profiles {:dev {:dependencies [[criterium/criterium "0.4.3"]
                                  [net.mikera/core.matrix.testing "0.0.4"]
                                  [expectations "2.1.1"]]}})
