(defproject clatrix/clatrix "0.1.0"
  :description "Because using matrices in Clojure needs to not suck."
  :url "https://github.com/tel/clatrix"
  :license {:name "MIT License"
            :url "http://www.opensource.org/licenses/MIT"}
  :resource-paths ["native"]
  :plugins [[lein-expectations "0.0.8"]]
  :dev-dependencies [[lein-expectations "0.0.8"]
                     [expectations "1.4.16"]]
  :dependencies [[org.clojure/clojure "1.4.0"]
                 [slingshot "0.8.0"]
                 [jblas/jblas "1.2.1"]
                 [jblas/native "1.2.0"]]
  :profiles {:dev {:dependencies [[expectations "1.4.16"]]}})
